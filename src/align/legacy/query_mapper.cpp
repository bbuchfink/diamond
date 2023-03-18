/****
DIAMOND protein aligner
Copyright (C) 2016-2021 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
						
Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#include <memory>
#include <algorithm>
#include "query_mapper.h"
#include "../../data/reference.h"
#include "../../dp/ungapped.h"
#include "../../output/output.h"
#include "../../output/output_format.h"
#include "../../output/daa/daa_write.h"
#include "../../output/target_culling.h"
#include "../../util/util.h"

using std::tie;
using std::list;
using std::vector;
using std::unique_ptr;
using std::min;
using std::max;
using std::set;
using std::string;

bool Target::envelopes(const ApproxHsp &t, double p) const
{
	for (list<ApproxHsp>::const_iterator i = ts.begin(); i != ts.end(); ++i)
		if (t.query_source_range.overlap_factor(i->query_source_range) >= p)
			return true;
	return false;
}

bool Target::is_enveloped(const Target &t, double p) const
{
	for (list<ApproxHsp>::const_iterator i = ts.begin(); i != ts.end(); ++i)
		if (!t.envelopes(*i, p))
			return false;
	return true;
}

bool Target::is_enveloped(PtrVector<Target>::const_iterator begin, PtrVector<Target>::const_iterator end, double p, int min_score) const
{
	for (; begin != end; ++begin)
		if (is_enveloped(**begin, p) && (*begin)->filter_score >= min_score)
			return true;
	return false;
}

void Target::add_ranges(vector<int32_t> &v) {
	for (const Hsp &hsp : hsps) {
		const int i0 = hsp.query_source_range.begin_ / INTERVAL,
			i1 = min(hsp.query_source_range.end_ / INTERVAL, int(v.size() - 1));
		for (int i = i0; i <= i1; ++i)
			v[i] = max(v[i], hsp.score);
	}
}

bool Target::is_outranked(const vector<int32_t> &v, double treshold) {
	for (const Hsp &hsp : hsps) {
		const int i0 = hsp.query_source_range.begin_ / INTERVAL,
			i1 = min(hsp.query_source_range.end_ / INTERVAL, int(v.size() - 1));
		for (int i = i0; i <= i1; ++i)
			if (hsp.score >= v[i] * treshold)
				return false;
	}
	return true;
}

QueryMapper::QueryMapper(size_t query_id, Search::Hit* begin, Search::Hit* end, const Search::Config &metadata) :
	source_hits(std::make_pair(begin, end)),
	query_id((unsigned)query_id),
	targets_finished(0),
	next_target(0),
	source_query_len(metadata.query->source_len((unsigned)query_id)),
	translated_query(metadata.query->translated(query_id)),
	metadata(metadata),
	target_parallel(false)
{
	seed_hits.reserve(source_hits.second - source_hits.first);
}

void QueryMapper::init()
{
	if(config.log_query)
		message_stream << "Query = " << metadata.query->ids()[query_id] << '\t' << query_id << std::endl;
	if (Stats::CBS::hauser(config.comp_based_stats))
		for (int i = 0; i < align_mode.query_contexts; ++i)
			query_cb.emplace_back(query_seq(i));
	targets.resize(count_targets());
	if (targets.empty())
		return;
	load_targets();
}

unsigned QueryMapper::count_targets()
{
	std::sort(source_hits.first, source_hits.second, Search::Hit::CmpSubject());
	const size_t n = source_hits.second - source_hits.first;
	const Search::Hit* hits = source_hits.first;
	size_t subject_id = std::numeric_limits<size_t>::max();
	unsigned n_subject = 0;
	for (size_t i = 0; i < n; ++i) {
		std::pair<size_t, size_t> l = metadata.target->seqs().local_position((uint64_t)hits[i].subject_);
		const unsigned frame = hits[i].query_ % align_mode.query_contexts;
		/*const Diagonal_segment d = config.comp_based_stats ? xdrop_ungapped(query_seq(frame), query_cb[frame], ref_seqs::get()[l.first], hits[i].seed_offset_, (int)l.second)
			: xdrop_ungapped(query_seq(frame), ref_seqs::get()[l.first], hits[i].seed_offset_, (int)l.second);*/
		if (target_parallel) {
			seed_hits.emplace_back(frame, (unsigned)l.first, (unsigned)l.second, (unsigned)hits[i].seed_offset_, DiagonalSegment());
			if (l.first != subject_id) {
				subject_id = l.first;
				++n_subject;
			}
		}
		else {
			const DiagonalSegment d = xdrop_ungapped(query_seq(frame), nullptr, metadata.target->seqs()[l.first], hits[i].seed_offset_, (int)l.second, false);
			if (d.score > 0) {
				if (l.first != subject_id) {
					subject_id = l.first;
					++n_subject;
				}
				seed_hits.emplace_back(frame, (unsigned)l.first, (unsigned)l.second, (unsigned)hits[i].seed_offset_, d);
			}
		}
	}
	return n_subject;
}

void QueryMapper::load_targets()
{
	unsigned subject_id = std::numeric_limits<unsigned>::max(), n = 0;
	for (size_t i = 0; i < seed_hits.size(); ++i) {
		if (seed_hits[i].subject_ != subject_id) {
			if (n > 0) {
				targets[n - 1].end = i;
			}
			const size_t oid = metadata.target->block_id2oid(seed_hits[i].subject_);
			targets.get(n) = new Target(i,
				seed_hits[i].subject_,
				metadata.target->seqs()[seed_hits[i].subject_],
				config.taxon_k ? metadata.db->taxon_nodes().rank_taxid(metadata.db->taxids(oid), Rank::species) : set<TaxId>());
			++n;
			subject_id = seed_hits[i].subject_;
		}
	}
	targets[n - 1].end = seed_hits.size();
}

void QueryMapper::rank_targets(double ratio, double factor, const int64_t max_target_seqs)
{
	if (config.taxon_k && config.toppercent == 100.0)
		return;
	std::stable_sort(targets.begin(), targets.end(), Target::compare_score);

	int score = 0;
	if (config.toppercent < 100) {
		score = int((double)targets[0].filter_score * (1.0 - config.toppercent / 100.0) * ratio);
	}
	else {
		int64_t min_idx = std::min((int64_t)targets.size(), max_target_seqs);
		score = int((double)targets[min_idx - 1].filter_score * ratio);
	}

	const int64_t cap = (config.toppercent < 100 || max_target_seqs == INT64_MAX) ? INT64_MAX : int64_t(max_target_seqs * factor);
	int64_t i = 0;
	for (; i < (int64_t)targets.size(); ++i)
		if (targets[i].filter_score < score || i >= cap)
			break;

	targets.erase(targets.begin() + i, targets.end());
}

void QueryMapper::score_only_culling(const int64_t max_target_seqs)
{
	static const double COV_INCLUDE_CUTOFF = 0.1;
	std::stable_sort(targets.begin(), targets.end(), config.toppercent == 100.0 ? Target::compare_evalue : Target::compare_score);
	unique_ptr<TargetCulling> target_culling(TargetCulling::get(max_target_seqs));
	const unsigned query_len = (unsigned)query_seq(0).length();
	PtrVector<Target>::iterator i;
	for (i = targets.begin(); i<targets.end();) {
		if (!score_matrix.report_cutoff((*i)->filter_score, (*i)->filter_evalue))
			break;
		int code;
		double cov;
		tie(code, cov) = target_culling->cull(**i);
		if (code == TargetCulling::FINISHED)
			break;
		else if (code == TargetCulling::NEXT) {
			i = targets.erase(i, i + 1);
		}
		else {
			if (cov < COV_INCLUDE_CUTOFF)
				target_culling->add(**i);
			++i;
		}
	}
	targets.erase(i, targets.end());
}

bool QueryMapper::generate_output(TextBuffer &buffer, Statistics &stat, const Search::Config& cfg)
{
	std::stable_sort(targets.begin(), targets.end(), config.toppercent == 100.0 ? Target::compare_evalue : Target::compare_score);

	unsigned n_hsp = 0, n_target_seq = 0, hit_hsps = 0;
	unique_ptr<TargetCulling> target_culling(TargetCulling::get(cfg.max_target_seqs));
	const unsigned query_len = (unsigned)query_seq(0).length();
	size_t seek_pos = 0;
	const char *query_title = metadata.query->ids()[query_id];
	unique_ptr<OutputFormat> f(cfg.output_format->clone());
	Output::Info info{ cfg.query->seq_info(query_id), true, cfg.db.get(), buffer, {} };

	for (size_t i = 0; i < targets.size(); ++i) {

		const BlockId subject_id = targets[i].subject_block_id;
		const OId database_id = metadata.target->block_id2oid(subject_id);
		string target_title;
		size_t dict_id;
		if (!cfg.blocked_processing)
			target_title = metadata.target->has_ids() ? metadata.target->ids()[subject_id] : metadata.db->seqid(database_id);
		else
			dict_id = metadata.target->dict_id(cfg.current_ref_block, subject_id, *metadata.db);
		const unsigned subject_len = (unsigned)metadata.target->seqs()[subject_id].length();
		targets[i].apply_filters(source_query_len, subject_len, query_title);
		if (targets[i].hsps.size() == 0)
			continue;

		const int c = target_culling->cull(targets[i]).first;
		if (c == TargetCulling::NEXT)
			continue;
		else if (c == TargetCulling::FINISHED)
			break;

		target_culling->add(targets[i]);
		
		hit_hsps = 0;
		for (list<Hsp>::iterator j = targets[i].hsps.begin(); j != targets[i].hsps.end(); ++j) {
			info.unaligned = false;
			if (config.max_hsps > 0 && hit_hsps >= config.max_hsps)
				break;

			if (cfg.blocked_processing) {
				if (n_hsp == 0)
					seek_pos = IntermediateRecord::write_query_intro(buffer, query_id);
				IntermediateRecord::write(buffer, *j, query_id, dict_id, database_id, cfg.output_format.get());
			}
			else {
				if (n_hsp == 0) {
					if (*f == OutputFormat::daa)
						seek_pos = write_daa_query_record(buffer, query_title, align_mode.query_translated ? metadata.query->source_seqs()[query_id] : metadata.query->seqs()[query_id]);
					else
						f->print_query_intro(info);
				}
				if (*f == OutputFormat::daa)
					write_daa_record(buffer, *j, safe_cast<uint32_t>(metadata.target->dict_id(cfg.current_ref_block, subject_id, *metadata.db)));
				else
					f->print_match(HspContext(*j,
						query_id,
						cfg.query->block_id2oid(query_id),
						translated_query,
						query_title,
						database_id,
						subject_len,
						target_title.c_str(),
						n_target_seq,
						hit_hsps,
						metadata.target->seqs()[subject_id]), info);
			}

			++n_hsp;
			++hit_hsps;
		}
		++n_target_seq;
	}

	if (n_hsp > 0) {
		if (!cfg.blocked_processing) {
			if (*f == OutputFormat::daa)
				finish_daa_query_record(buffer, seek_pos);
			else
				f->print_query_epilog(info);
		}
		else
			IntermediateRecord::finish_query(buffer, seek_pos);
	}
	else if (!cfg.blocked_processing && *f != OutputFormat::daa && config.report_unaligned != 0) {
		f->print_query_intro(info);
		f->print_query_epilog(info);
	}

	if (!cfg.blocked_processing) {
		stat.inc(Statistics::MATCHES, n_hsp);
		stat.inc(Statistics::PAIRWISE, n_target_seq);
		if (n_hsp > 0)
			stat.inc(Statistics::ALIGNED);
	}
	
	return n_hsp > 0;
}

void Target::inner_culling()
{
	hsps.sort();
	if (hsps.size() > 0) {
		filter_score = hsps.front().score;
		filter_evalue = hsps.front().evalue;
	}
	else {
		filter_score = 0;
		filter_evalue = DBL_MAX;
	}
	for (list<Hsp>::iterator i = hsps.begin(); i != hsps.end();) {
		if (i->query_range_enveloped_by(hsps.begin(), i, 0.5))
			i = hsps.erase(i);
		else
			++i;
	}
}

void Target::apply_filters(int dna_len, int subject_len, const char *query_title)
{
	for (list<Hsp>::iterator i = hsps.begin(); i != hsps.end();) {
		if (i->id_percent() < config.min_id
			|| i->query_cover_percent(dna_len) < config.query_cover
			|| i->subject_cover_percent(subject_len) < config.subject_cover)
			i = hsps.erase(i);
		else
			++i;
	}
}
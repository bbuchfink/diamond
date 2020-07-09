/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

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
#include "../extend_ungapped.h"
#include "../../output/output.h"
#include "../../output/output_format.h"
#include "../../output/daa_write.h"
#include "../../output/target_culling.h"

using namespace std;

bool Target::envelopes(const Hsp_traits &t, double p) const
{
	for (list<Hsp_traits>::const_iterator i = ts.begin(); i != ts.end(); ++i)
		if (t.query_source_range.overlap_factor(i->query_source_range) >= p)
			return true;
	return false;
}

bool Target::is_enveloped(const Target &t, double p) const
{
	for (list<Hsp_traits>::const_iterator i = ts.begin(); i != ts.end(); ++i)
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

int QueryMapper::raw_score_cutoff() const
{
	return score_matrix.rawscore(config.min_bit_score == 0 ? score_matrix.bitscore(config.max_evalue, (unsigned)query_seq(0).length()) : config.min_bit_score);
}

QueryMapper::QueryMapper(const Parameters &params, size_t query_id, hit* begin, hit* end, const Metadata &metadata, bool target_parallel) :
	parameters(params),
	source_hits(std::make_pair(begin, end)),
	query_id((unsigned)query_id),
	targets_finished(0),
	next_target(0),
	source_query_len(get_source_query_len((unsigned)query_id)),
	translated_query(get_translated_query(query_id)),
	target_parallel(target_parallel),
	metadata(metadata)
{
	seed_hits.reserve(source_hits.second - source_hits.first);
}

void QueryMapper::init()
{
	if(config.log_query)
		cout << "Query = " << query_ids::get()[query_id] << '\t' << query_id << endl;
	if (config.comp_based_stats == 1)
		for (unsigned i = 0; i < align_mode.query_contexts; ++i)
			query_cb.emplace_back(query_seq(i));
	targets.resize(count_targets());
	if (targets.empty())
		return;
	load_targets();
}

unsigned QueryMapper::count_targets()
{
	std::sort(source_hits.first, source_hits.second, hit::CmpSubject());
	const size_t n = source_hits.second - source_hits.first;
	const hit* hits = source_hits.first;
	size_t subject_id = std::numeric_limits<size_t>::max();
	unsigned n_subject = 0;
	for (size_t i = 0; i < n; ++i) {
		std::pair<size_t, size_t> l = ref_seqs::data_->local_position((uint64_t)hits[i].subject_);
		const unsigned frame = hits[i].query_ % align_mode.query_contexts;
		/*const Diagonal_segment d = config.comp_based_stats ? xdrop_ungapped(query_seq(frame), query_cb[frame], ref_seqs::get()[l.first], hits[i].seed_offset_, (int)l.second)
			: xdrop_ungapped(query_seq(frame), ref_seqs::get()[l.first], hits[i].seed_offset_, (int)l.second);*/
		if (target_parallel) {
			seed_hits.emplace_back(frame, (unsigned)l.first, (unsigned)l.second, (unsigned)hits[i].seed_offset_, Diagonal_segment());
			if (l.first != subject_id) {
				subject_id = l.first;
				++n_subject;
			}
		}
		else {
			const Diagonal_segment d = xdrop_ungapped(query_seq(frame), ref_seqs::get()[l.first], hits[i].seed_offset_, (int)l.second);
			if (d.score >= config.min_ungapped_raw_score) {
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
			targets.get(n) = new Target(i,
				seed_hits[i].subject_,
				config.taxon_k ? metadata.taxon_nodes->rank_taxid((*metadata.taxon_list)[ReferenceDictionary::get().block_to_database_id(seed_hits[i].subject_)], Rank::species) : set<unsigned>());
			++n;
			subject_id = seed_hits[i].subject_;
		}
	}
	targets[n - 1].end = seed_hits.size();
}

void QueryMapper::rank_targets(double ratio, double factor)
{
	if (config.taxon_k && config.toppercent == 100.0)
		return;
	std::stable_sort(targets.begin(), targets.end(), Target::compare);

	int score = 0;
	if (config.toppercent < 100) {
		score = int((double)targets[0].filter_score * (1.0 - config.toppercent / 100.0) * ratio);
	}
	else {
		size_t min_idx = std::min(targets.size(), (size_t)config.max_alignments);
		score = int((double)targets[min_idx - 1].filter_score * ratio);
	}

	const size_t cap = (config.toppercent < 100 || config.max_alignments == std::numeric_limits<uint64_t>::max()) ? std::numeric_limits<uint64_t>::max() : size_t(config.max_alignments*factor);
	size_t i = 0;
	for (; i < targets.size(); ++i)
		if (targets[i].filter_score < score || i >= cap)
			break;

	if (config.benchmark_ranking)
		for (size_t j = i; j < targets.size(); ++j)
			targets[j].outranked = true;
	else
		targets.erase(targets.begin() + i, targets.end());
}

void QueryMapper::score_only_culling()
{
	std::stable_sort(targets.begin(), targets.end(), Target::compare);
	unique_ptr<TargetCulling> target_culling(TargetCulling::get());
	const unsigned query_len = (unsigned)query_seq(0).length();
	PtrVector<Target>::iterator i;
	for (i = targets.begin(); i<targets.end();) {
		if ((config.min_bit_score == 0 && score_matrix.evalue((*i)->filter_score, query_len) > config.max_evalue)
			|| score_matrix.bitscore((*i)->filter_score) < config.min_bit_score)
			break;
		const int c = target_culling->cull(**i);
		if (c == TargetCulling::FINISHED)
			break;
		else if (c == TargetCulling::NEXT) {
			if (config.benchmark_ranking)
				(*i++)->outranked = true;
			else
				i = targets.erase(i, i + 1);
		}
		else {
			target_culling->add(**i);
			++i;
		}
	}
	targets.erase(i, targets.end());
}

bool QueryMapper::generate_output(TextBuffer &buffer, Statistics &stat)
{
	std::stable_sort(targets.begin(), targets.end(), Target::compare);

	unsigned n_hsp = 0, n_target_seq = 0, hit_hsps = 0;
	unique_ptr<TargetCulling> target_culling(TargetCulling::get());
	const unsigned query_len = (unsigned)query_seq(0).length();
	size_t seek_pos = 0;
	const char *query_title = query_ids::get()[query_id];
	unique_ptr<Output_format> f(output_format->clone());

	for (size_t i = 0; i < targets.size(); ++i) {

		if ((config.min_bit_score == 0 && score_matrix.evalue(targets[i].filter_score, query_len) > config.max_evalue)
			|| score_matrix.bitscore(targets[i].filter_score) < config.min_bit_score)
			break;

		const size_t subject_id = targets[i].subject_block_id;
		const unsigned database_id = ReferenceDictionary::get().block_to_database_id(subject_id);
		const unsigned subject_len = (unsigned)ref_seqs::get()[subject_id].length();
		const char *ref_title = ref_ids::get()[subject_id];
		targets[i].apply_filters(source_query_len, subject_len, query_title, ref_title);
		if (targets[i].hsps.size() == 0)
			continue;

		const int c = target_culling->cull(targets[i]);
		if (c == TargetCulling::NEXT)
			continue;
		else if (c == TargetCulling::FINISHED)
			break;

		if (targets[i].outranked)
			stat.inc(Statistics::OUTRANKED_HITS);
		target_culling->add(targets[i]);
		
		hit_hsps = 0;
		for (list<Hsp>::iterator j = targets[i].hsps.begin(); j != targets[i].hsps.end(); ++j) {
			if (config.max_hsps > 0 && hit_hsps >= config.max_hsps)
				break;

			if (blocked_processing) {
				if (n_hsp == 0)
					seek_pos = IntermediateRecord::write_query_intro(buffer, query_id);
				IntermediateRecord::write(buffer, *j, query_id, subject_id);
			}
			else {
				if (n_hsp == 0) {
					if (*f == Output_format::daa)
						seek_pos = write_daa_query_record(buffer, query_title, align_mode.query_translated ? query_source_seqs::get()[query_id] : query_seqs::get()[query_id]);
					else
						f->print_query_intro(query_id, query_title, source_query_len, buffer, false);
				}
				if (*f == Output_format::daa)
					write_daa_record(buffer, *j, subject_id);
				else
					f->print_match(Hsp_context(*j,
						query_id,
						translated_query,
						query_title,
						subject_id,
						database_id,
						ref_title,
						subject_len,
						n_target_seq,
						hit_hsps,
						ref_seqs::get()[subject_id]), metadata, buffer);
			}

			++n_hsp;
			++hit_hsps;
		}
		++n_target_seq;
	}

	if (n_hsp > 0) {
		if (!blocked_processing) {
			if (*f == Output_format::daa)
				finish_daa_query_record(buffer, seek_pos);
			else
				f->print_query_epilog(buffer, query_title, false, parameters);
		}
		else
			IntermediateRecord::finish_query(buffer, seek_pos);
	}
	else if (!blocked_processing && *f != Output_format::daa && config.report_unaligned != 0) {
		f->print_query_intro(query_id, query_title, source_query_len, buffer, true);
		f->print_query_epilog(buffer, query_title, true, parameters);
	}

	if (!blocked_processing) {
		stat.inc(Statistics::MATCHES, n_hsp);
		stat.inc(Statistics::PAIRWISE, n_target_seq);
		if (n_hsp > 0)
			stat.inc(Statistics::ALIGNED);
	}
	
	return n_hsp > 0;
}

void Target::inner_culling(int cutoff)
{
	hsps.sort();
	if (hsps.size() > 0)
		filter_score = hsps.front().score;
	else
		filter_score = 0;
	for (list<Hsp>::iterator i = hsps.begin(); i != hsps.end();) {
		if (i->is_enveloped_by(hsps.begin(), i, 0.5) || (int)i->score < cutoff)
			i = hsps.erase(i);
		else
			++i;
	}
}

void Target::apply_filters(int dna_len, int subject_len, const char *query_title, const char *ref_title)
{
	for (list<Hsp>::iterator i = hsps.begin(); i != hsps.end();) {
		if (i->id_percent() < config.min_id
			|| i->query_cover_percent(dna_len) < config.query_cover
			|| i->subject_cover_percent(subject_len) < config.subject_cover
			|| (config.no_self_hits
				&& i->identities == i->length
				&& i->query_source_range.length() == (int)dna_len
				&& i->subject_range.length() == (int)subject_len
				&& strcmp(query_title, ref_title) == 0)
			|| (config.filter_locus && !i->subject_range.includes(config.filter_locus)))
			i = hsps.erase(i);
		else
			++i;
	}
}
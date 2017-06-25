/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#include "align_queries.h"
#include "query_mapper.h"
#include "../data/reference.h"
#include "extend_ungapped.h"
#include "../output/output.h"
#include "../output/output_format.h"
#include "../output/daa_write.h"

Query_mapper::Query_mapper() :
	source_hits(get_query_data()),
	query_id(source_hits.first->query_ / align_mode.query_contexts),
	targets_finished(0),
	next_target(0),
	source_query_len(get_source_query_len(query_id)),
	unaligned_from(query_queue.last_query+1)
{	
	query_queue.last_query = query_id;
	seed_hits.reserve(source_hits.second - source_hits.first);
}

int Query_mapper::raw_score_cutoff() const
{
	return score_matrix.rawscore(config.min_bit_score == 0 ? score_matrix.bitscore(config.max_evalue, ref_header.letters, (unsigned)query_seq(0).length()) : config.min_bit_score);
}

Query_mapper::Query_mapper(size_t query_id, Trace_pt_list::iterator begin, Trace_pt_list::iterator end) :
	source_hits(std::make_pair(begin, end)),
	query_id((unsigned)query_id),
	targets_finished(0),
	next_target(0),
	source_query_len(get_source_query_len((unsigned)query_id))
{
	seed_hits.reserve(source_hits.second - source_hits.first);
}

void Query_mapper::init()
{
	if(config.log_query)
		cout << "Query = " << query_ids::get()[query_id].c_str() << endl;
	if (config.comp_based_stats == 1)
		for (unsigned i = 0; i < align_mode.query_contexts; ++i)
			query_cb.push_back(Bias_correction(query_seq(i)));
	if (config.ext == Config::greedy || config.ext == Config::more_greedy)
		for (unsigned i = 0; i < align_mode.query_contexts; ++i)
			profile.push_back(Long_score_profile(query_seq(i)));
			//profile.push_back(Long_score_profile());
	targets.resize(count_targets());
	if (targets.empty())
		return;
	load_targets();
	if (config.ext == Config::floating_xdrop)
		rank_targets(config.rank_ratio == -1 ? 0.6 : config.rank_ratio);
}

pair<Trace_pt_list::iterator, Trace_pt_list::iterator> Query_mapper::get_query_data()
{
	const Trace_pt_list::iterator begin = query_queue.trace_pt_pos;
	if (begin == query_queue.trace_pt_end)
		return pair<Trace_pt_list::iterator, Trace_pt_list::iterator>(begin, begin);
	const unsigned c = align_mode.query_contexts, query = begin->query_ / c;
	Trace_pt_list::iterator end = begin;
	for (; end < query_queue.trace_pt_end && end->query_ / c == query; ++end);
	query_queue.trace_pt_pos = end;
	return pair<Trace_pt_list::iterator, Trace_pt_list::iterator>(begin, end);
}

unsigned Query_mapper::count_targets()
{
	std::sort(source_hits.first, source_hits.second, hit::cmp_subject);
	const size_t n = source_hits.second - source_hits.first;
	const Trace_pt_list::iterator hits = source_hits.first;
	size_t subject_id = std::numeric_limits<size_t>::max();
	unsigned n_subject = 0;
	for (size_t i = 0; i < n; ++i) {
		std::pair<size_t, size_t> l = ref_seqs::data_->local_position(hits[i].subject_);
		const unsigned frame = hits[i].query_ % align_mode.query_contexts;
		/*const Diagonal_segment d = config.comp_based_stats ? xdrop_ungapped(query_seq(frame), query_cb[frame], ref_seqs::get()[l.first], hits[i].seed_offset_, (int)l.second)
			: xdrop_ungapped(query_seq(frame), ref_seqs::get()[l.first], hits[i].seed_offset_, (int)l.second);*/
		const Diagonal_segment d = xdrop_ungapped(query_seq(frame), ref_seqs::get()[l.first], hits[i].seed_offset_, (int)l.second);
		if (d.score >= config.min_ungapped_raw_score) {
			if (l.first != subject_id) {
				subject_id = l.first;
				++n_subject;
			}
			seed_hits.push_back(Seed_hit(frame,
				(unsigned)l.first,
				(unsigned)l.second,
				hits[i].seed_offset_,
				d));
		}
	}
	return n_subject;
}

void Query_mapper::load_targets()
{
	unsigned subject_id = std::numeric_limits<unsigned>::max(), n = 0;
	for (size_t i = 0; i < seed_hits.size(); ++i) {
		if (seed_hits[i].subject_ != subject_id) {
			if (n > 0) {
				targets[n - 1].end = i;
				get_prefilter_score(n - 1);
			}
			targets.get(n) = new Target(i, seed_hits[i].subject_);
			++n;
			subject_id = seed_hits[i].subject_;
		}
	}
	targets[n - 1].end = seed_hits.size();
	get_prefilter_score(n - 1);
}

void Query_mapper::rank_targets(double ratio)
{
	std::sort(targets.begin(), targets.end(), Target::compare);

	int score = 0;
	if (config.toppercent < 100) {
		score = int((double)targets[0].filter_score * (1.0 - config.toppercent / 100.0) * ratio);
	}
	else {
		size_t min_idx = std::min(targets.size(), (size_t)config.max_alignments);
		score = int((double)targets[min_idx - 1].filter_score * ratio);
	}

	unsigned i = 0;
	for (; i < targets.size(); ++i)
		if (targets[i].filter_score < score)
			break;
	
	if (config.benchmark_ranking)
		for (unsigned j = i; j < targets.size(); ++j)
			targets[j].outranked = true;
	else
		targets.erase(targets.begin() + i, targets.end());
}

bool Query_mapper::generate_output(Text_buffer &buffer, Statistics &stat)
{
	std::sort(targets.begin(), targets.end(), Target::compare);

	unsigned n_hsp = 0, n_target_seq = 0, hit_hsps = 0;
	const unsigned top_score = targets.empty() ? 0 : targets[0].filter_score, query_len = (unsigned)query_seq(0).length();
	size_t seek_pos = 0;
	const char *query_title = query_ids::get()[query_id].c_str();
	auto_ptr<Output_format> f(output_format->clone());

	for (size_t i = 0; i < targets.size(); ++i) {
		if ((config.min_bit_score == 0 && score_matrix.evalue(targets[i].filter_score, config.db_size, query_len) > config.max_evalue)
			|| score_matrix.bitscore(targets[i].filter_score) < config.min_bit_score)
			break;

		if (!config.output_range(n_target_seq, targets[i].filter_score, top_score))
			break;

		if (targets[i].outranked)
			stat.inc(Statistics::OUTRANKED_HITS);

		const unsigned subject_len = (unsigned)ref_seqs::get()[targets[i].subject_id].length();
		
		hit_hsps = 0;
		for (list<Hsp_data>::iterator j = targets[i].hsps.begin(); j != targets[i].hsps.end(); ++j) {
			if (hit_hsps >= config.max_hsps)
				break;
			const char *ref_title = ref_ids::get()[targets[i].subject_id].c_str();
			if (j->id_percent() < config.min_id
				|| j->query_cover_percent(source_query_len) < config.query_cover
				|| j->subject_cover_percent(subject_len) < config.subject_cover
				|| (config.no_self_hits &&
					(config.ext == Config::more_greedy || (j->identities == j->length && j->query_source_range.length() == (int)source_query_len && j->subject_range.length() == (int)subject_len))
					&& strcmp(query_title, ref_title) == 0))
				continue;

			if (blocked_processing) {
				if (n_hsp == 0)
					seek_pos = Intermediate_record::write_query_intro(buffer, query_id);
				Intermediate_record::write(buffer, *j, query_id, targets[i].subject_id);
			}
			else {
				if (n_hsp == 0) {
					if (*f == Output_format::daa)
						seek_pos = write_daa_query_record(buffer, query_title, align_mode.query_translated ? query_source_seqs::get()[query_id] : query_seqs::get()[query_id]);
					else
						f->print_query_intro(query_id, query_title, source_query_len, buffer, false);
				}
				if (*f == Output_format::daa)
					write_daa_record(buffer, *j, query_id, targets[i].subject_id);
				else
					f->print_match(Hsp_context(*j,
						query_id,
						query_seq(j->frame),
						query_source_seq(),
						query_title,
						targets[i].subject_id,
						targets[i].subject_id,
						ref_title,
						subject_len,
						n_target_seq,
						hit_hsps), buffer);
			}

			if(hit_hsps == 0)
				++n_target_seq;
			++n_hsp;
			++hit_hsps;
			if (config.alignment_traceback && j->gap_openings > 0)
				stat.inc(Statistics::GAPPED);
			stat.inc(Statistics::SCORE_TOTAL, j->score);
		}
	}

	if (n_hsp > 0) {
		if (!blocked_processing) {
			if (*f == Output_format::daa)
				finish_daa_query_record(buffer, seek_pos);
			else
				f->print_query_epilog(buffer, query_title, false);
		}
		else
			Intermediate_record::finish_query(buffer, seek_pos);
	}
	else if (!blocked_processing && *f != Output_format::daa && config.report_unaligned != 0) {
		f->print_query_intro(query_id, query_title, source_query_len, buffer, true);
		f->print_query_epilog(buffer, query_title, true);
	}

	if (!blocked_processing) {
		stat.inc(Statistics::MATCHES, n_hsp);
		stat.inc(Statistics::PAIRWISE, n_target_seq);
		if (n_hsp > 0)
			stat.inc(Statistics::ALIGNED);
	}
	
	return n_hsp > 0;
}
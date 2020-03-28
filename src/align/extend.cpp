/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen

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

#include "extend.h"
#include "../data/queries.h"
#include "../basic/config.h"
#include "../dp/comp_based_stats.h"
#include "target.h"
#include "../dp/dp.h"
#include "../util/log_stream.h"
#include "../data/reference.h"

using std::vector;
using std::list;
using std::array;

namespace Extension {

void load_hits(Trace_pt_list::iterator begin, Trace_pt_list::iterator end, FlatArray<SeedHit> &hits, vector<size_t> &target_block_ids) {
	hits.clear();
	target_block_ids.clear();
	if (begin >= end)
		return;
	std::sort(begin, end, hit::cmp_subject);
	size_t target = SIZE_MAX;
	for (Trace_pt_list::iterator i = begin; i < end; ++i) {
		std::pair<size_t, size_t> l = ref_seqs::data_->local_position(i->subject_);
		if (l.first != target) {
			hits.next();
			target = l.first;
			target_block_ids.push_back(target);
		}
		hits.push_back({ (int)i->seed_offset_, (int)l.second, i->query_ % align_mode.query_contexts });
	}
}

vector<Match> extend(const Parameters &params, size_t query_id, Trace_pt_list::iterator begin, Trace_pt_list::iterator end, const Metadata &metadata, Statistics &stat, int flags) {
	const unsigned contexts = align_mode.query_contexts;
	vector<sequence> query_seq;
	vector<Bias_correction> query_cb;

	if (config.log_query || flags & TARGET_PARALLEL)
		log_stream << "Query = " << query_ids::get()[query_id] << endl;

	for (unsigned i = 0; i < contexts; ++i)
		query_seq.push_back(query_seqs::get()[query_id*contexts + i]);

	task_timer timer(flags & TARGET_PARALLEL ? 3 : UINT_MAX);
	if (config.comp_based_stats == 1) {
		timer.go("Computing CBS");
		for (unsigned i = 0; i < contexts; ++i)
			query_cb.emplace_back(query_seq[i]);
		timer.finish();
	}

	const int source_query_len = align_mode.query_translated ? (int)query_source_seqs::get()[query_id].length() : (int)query_seqs::get()[query_id].length();

	timer.go("Loading seed hits");
	thread_local FlatArray<SeedHit> seed_hits;
	thread_local vector<size_t> target_block_ids;
	load_hits(begin, end, seed_hits, target_block_ids);
	stat.inc(Statistics::TARGET_HITS0, target_block_ids.size());

	if (config.gapped_filter_score > 0.0 || config.gapped_filter_evalue > 0.0) {
		timer.go("Computing gapped filter");
		gapped_filter(query_seq.data(), query_cb.data(), seed_hits, target_block_ids, stat);
	}
	stat.inc(Statistics::TARGET_HITS1, target_block_ids.size());

	timer.go("Computing chaining");
	vector<WorkTarget> targets = ungapped_stage(query_seq.data(), query_cb.data(), seed_hits, target_block_ids.data(), flags);
	stat.inc(Statistics::TARGET_HITS2, targets.size());
	
	timer.go("Computing ranking");
	rank_targets(targets, config.rank_ratio == -1 ? (query_seq[0].length() > 50 ? 0.6 : 0.9) : config.rank_ratio, config.rank_factor == -1 ? 1e3 : config.rank_factor);
	stat.inc(Statistics::TARGET_HITS3, targets.size());
	timer.finish();

	/*if (config.gapped_filter_score > 0.0 || config.gapped_filter_evalue > 0.0)
		targets = gapped_filter(query_seq.data(), query_cb.data(), targets, stat);
	stat.inc(Statistics::TARGET_HITS2, targets.size());*/

	vector<Target> aligned_targets = align(targets, query_seq.data(), query_cb.data(), flags, stat);
	timer.go("Computing score only culling");
	score_only_culling(aligned_targets);
	stat.inc(Statistics::TARGET_HITS4, aligned_targets.size());
	timer.finish();

	vector<Match> matches = align(aligned_targets, query_seq.data(), query_cb.data(), source_query_len, flags, stat);
	timer.go("Computing culling");
	culling(matches, source_query_len, query_ids::get()[query_id]);
	stat.inc(Statistics::TARGET_HITS5, matches.size());

	return matches;
}

}
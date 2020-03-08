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

using std::vector;
using std::list;
using std::array;

namespace Extension {

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

	vector<WorkTarget> targets = ungapped_stage(query_seq.data(), query_cb.data(), begin, end, flags);
	stat.inc(Statistics::TARGET_HITS0, targets.size());

	timer.go("Computing ranking");
	rank_targets(targets, config.rank_ratio == -1 ? (query_seq[0].length() > 50 ? 0.6 : 0.9) : config.rank_ratio, config.rank_factor == -1 ? 1e3 : config.rank_factor);
	stat.inc(Statistics::TARGET_HITS1, targets.size());
	timer.finish();

	vector<Target> aligned_targets = align(targets, query_seq.data(), query_cb.data(), flags, stat);
	timer.go("Computing score only culling");
	score_only_culling(aligned_targets);
	stat.inc(Statistics::TARGET_HITS2, aligned_targets.size());
	timer.finish();

	vector<Match> matches = align(aligned_targets, query_seq.data(), query_cb.data(), source_query_len, flags, stat);
	timer.go("Computing culling");
	culling(matches, source_query_len, query_ids::get()[query_id]);
	stat.inc(Statistics::TARGET_HITS3, matches.size());

	return matches;
}

}
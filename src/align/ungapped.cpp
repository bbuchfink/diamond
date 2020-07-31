/****
DIAMOND protein aligner
Copyright (C) 2020 Max Planck Society for the Advancement of Science e.V.

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

#include <array>
#include <algorithm>
#include <limits.h>
#include <utility>
#include <atomic>
#include <mutex>
#include <limits.h>
#include "../basic/config.h"
#include "../dp/comp_based_stats.h"
#include "target.h"
#include "../dp/ungapped.h"
#include "target.h"
#include "../data/reference.h"
#include "../util/log_stream.h"
#include "../util/parallel/thread_pool.h"
#include "../chaining/chaining.h"

using std::array;
using std::vector;
using std::list;
using std::atomic;
using std::mutex;

namespace Extension {

WorkTarget ungapped_stage(SeedHit *begin, SeedHit *end, const sequence *query_seq, const Bias_correction *query_cb, uint32_t block_id) {
	array<vector<Diagonal_segment>, MAX_CONTEXT> diagonal_segments;
	WorkTarget target(block_id, ref_seqs::get()[block_id]);
	if (config.ext == "full" && (int)query_seq[0].length() < config.full_sw_len)
		return target;
	std::sort(begin, end);
	for (const SeedHit *hit = begin; hit < end; ++hit) {
		if (!diagonal_segments[hit->frame].empty() && diagonal_segments[hit->frame].back().diag() == hit->diag() && diagonal_segments[hit->frame].back().subject_end() >= hit->j)
			continue;
		const auto d = xdrop_ungapped(query_seq[hit->frame], target.seq, hit->i, hit->j);
		if (d.score > 0) {
			target.ungapped_score = std::max(target.ungapped_score, d.score);
			diagonal_segments[hit->frame].push_back(d);
		}
	}
	for (unsigned frame = 0; frame < align_mode.query_contexts; ++frame) {
		if (diagonal_segments[frame].empty())
			continue;
		std::stable_sort(diagonal_segments[frame].begin(), diagonal_segments[frame].end(), Diagonal_segment::cmp_diag);
		pair<int, list<Hsp_traits>> hsp = greedy_align(query_seq[frame], target.seq, diagonal_segments[frame].begin(), diagonal_segments[frame].end(), config.log_extend, frame);
		target.filter_score = std::max(target.filter_score, hsp.first);
		target.hsp[frame] = std::move(hsp.second);
		target.hsp[frame].sort(Hsp_traits::cmp_diag);
	}
	return target;
}

void ungapped_stage_worker(size_t i, size_t thread_id, const sequence *query_seq, const Bias_correction *query_cb, FlatArray<SeedHit> *seed_hits, const uint32_t*target_block_ids, vector<WorkTarget> *out, mutex *mtx) {
	WorkTarget target = ungapped_stage(seed_hits->begin(i), seed_hits->end(i), query_seq, query_cb, target_block_ids[i]);
	{
		std::lock_guard<mutex> guard(*mtx);
		out->push_back(std::move(target));
	}
}

vector<WorkTarget> ungapped_stage(const sequence *query_seq, const Bias_correction *query_cb, FlatArray<SeedHit> &seed_hits, const uint32_t *target_block_ids, int flags) {
	vector<WorkTarget> targets;
	if (seed_hits.size() == 0)
		return targets;
	if (flags & TARGET_PARALLEL) {
		mutex mtx;
		Util::Parallel::scheduled_thread_pool_auto(config.threads_, seed_hits.size(), ungapped_stage_worker, query_seq, query_cb, &seed_hits, target_block_ids, &targets, &mtx);
	}
	else
		for (size_t i = 0; i < seed_hits.size(); ++i)
			targets.push_back(ungapped_stage(seed_hits.begin(i), seed_hits.end(i), query_seq, query_cb, target_block_ids[i]));

	return targets;
}

}
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
#include "../util/data_structures/flat_array.h"
#include "../util/log_stream.h"
#include "../util/parallel/thread_pool.h"
#include "../chaining/chaining.h"

using std::array;
using std::vector;
using std::list;
using std::atomic;
using std::mutex;

namespace Extension {

struct SeedHit {
	int i, j;
	unsigned frame;
};

WorkTarget ungapped_stage(const SeedHit *begin, const SeedHit *end, const sequence *query_seq, const Bias_correction *query_cb, size_t block_id) {
	array<vector<Diagonal_segment>, MAX_CONTEXT> diagonal_segments;
	WorkTarget target(block_id, ref_seqs::get()[block_id]);
	for (const SeedHit *hit = begin; hit < end; ++hit) {
		const auto d = xdrop_ungapped(query_seq[hit->frame], target.seq, hit->i, hit->j);
		if(d.score >= config.min_ungapped_raw_score) diagonal_segments[hit->frame].push_back(d);
	}
	for (unsigned frame = 0; frame < align_mode.query_contexts; ++frame) {
		if (diagonal_segments[frame].empty())
			continue;
		std::stable_sort(diagonal_segments[frame].begin(), diagonal_segments[frame].end(), Diagonal_segment::cmp_diag);
		pair<int, list<Hsp_traits>> hsp = greedy_align(query_seq[frame], query_cb[frame], target.seq, diagonal_segments[frame].begin(), diagonal_segments[frame].end(), config.log_extend, frame);
		target.filter_score = std::max(target.filter_score, hsp.first);
		target.hsp[frame] = std::move(hsp.second);
	}
	return target;
}

void ungapped_stage_worker(size_t i, size_t thread_id, const sequence *query_seq, const Bias_correction *query_cb, const FlatArray<SeedHit> *seed_hits, size_t *target_block_ids, vector<WorkTarget> *out, mutex *mtx) {
	WorkTarget target = ungapped_stage(seed_hits->begin(i), seed_hits->end(i), query_seq, query_cb, target_block_ids[i]);
	{
		std::lock_guard<mutex> guard(*mtx);
		out->push_back(std::move(target));
	}
}

vector<WorkTarget> ungapped_stage(const sequence *query_seq, const Bias_correction *query_cb, Trace_pt_list::iterator begin, Trace_pt_list::iterator end, int flags) {
	task_timer timer("Loading seed hits", flags & TARGET_PARALLEL ? 3 : UINT_MAX);
	vector<WorkTarget> targets;
	if (begin >= end)
		return targets;
	std::sort(begin, end, hit::cmp_subject);
	size_t target = SIZE_MAX;
	thread_local FlatArray<SeedHit> hits;
	thread_local vector<size_t> target_block_ids;
	hits.clear();
	target_block_ids.clear();
	for (Trace_pt_list::iterator i = begin; i < end; ++i) {
		std::pair<size_t, size_t> l = ref_seqs::data_->local_position(i->subject_);
		if (l.first != target) {
			hits.next();
			target = l.first;
			target_block_ids.push_back(target);
		}
		hits.push_back({ (int)i->seed_offset_, (int)l.second, i->query_ % align_mode.query_contexts });
	}

	timer.go("Computing chaining");
	if (flags & TARGET_PARALLEL) {
		mutex mtx;
		Util::Parallel::scheduled_thread_pool_auto(config.threads_, hits.size(), ungapped_stage_worker, query_seq, query_cb, &hits, target_block_ids.data(), &targets, &mtx);
	}
	else
		for (size_t i = 0; i < hits.size(); ++i)
			targets.push_back(ungapped_stage(hits.begin(i), hits.end(i), query_seq, query_cb, target_block_ids[i]));

	return targets;
}

}
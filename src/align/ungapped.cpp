/****
DIAMOND protein aligner
Copyright (C) 2021 Max Planck Society for the Advancement of Science e.V.

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
#include "../stats/hauser_correction.h"
#include "target.h"
#include "../dp/ungapped.h"
#include "target.h"
#include "../data/reference.h"
#include "../util/log_stream.h"
#include "../util/parallel/thread_pool.h"
#include "../chaining/chaining.h"
#include "../dp/dp.h"
#include "def.h"
#include "../util/geo/geo.h"

using std::array;
using std::vector;
using std::list;
using std::atomic;
using std::mutex;
using std::pair;

namespace Extension {

std::vector<int16_t*> target_matrices;
std::mutex target_matrices_lock;
atomic<size_t> target_matrix_count(0);

WorkTarget::WorkTarget(BlockId block_id, const Sequence& seq, int query_len, const ::Stats::Composition& query_comp, const int16_t** query_matrix) :
	block_id(block_id),
	seq(seq),
	done(false)
{
	ungapped_score.fill(0);
	matrix = ::Stats::TargetMatrix(query_comp, query_len, seq);
}

WorkTarget ungapped_stage(FlatArray<SeedHit>::DataIterator begin, FlatArray<SeedHit>::DataIterator end, const Sequence *query_seq, const Bias_correction *query_cb, const ::Stats::Composition& query_comp, const int16_t** query_matrix, uint32_t block_id, Statistics& stat, const Block& targets, const Mode mode) {
	array<vector<DiagonalSegment>, MAX_CONTEXT> diagonal_segments;
	TaskTimer timer;
	const SequenceSet& ref_seqs = targets.seqs(), &ref_seqs_unmasked = targets.unmasked_seqs();
	const bool masking = config.comp_based_stats == ::Stats::CBS::COMP_BASED_STATS_AND_MATRIX_ADJUST ? ::Stats::use_seg_masking(query_seq[0], ref_seqs_unmasked[block_id]) : true;
	WorkTarget target(block_id, masking ? ref_seqs[block_id] : ref_seqs_unmasked[block_id], ::Stats::count_true_aa(query_seq[0]), query_comp, query_matrix);
	stat.inc(Statistics::TIME_MATRIX_ADJUST, timer.microseconds());
	if (target.adjusted_matrix())
		stat.inc(Statistics::MATRIX_ADJUST_COUNT);

	if (mode == Mode::FULL) {
		for (FlatArray<SeedHit>::DataIterator hit = begin; hit < end; ++hit)
			target.ungapped_score[hit->frame] = std::max(target.ungapped_score[hit->frame], hit->score);
		return target;
	}
	if (end - begin == 1 && align_mode.query_translated) { // ???
		target.ungapped_score[begin->frame] = begin->score;
		target.hsp[begin->frame].emplace_back(begin->diag(), begin->diag(), begin->score, begin->frame, begin->query_range(), begin->target_range(), begin->diag_segment());
		return target;
	}
	std::sort(begin, end);
	const bool with_diag_filter = config.hamming_ext || config.diag_filter_cov > 0 || config.diag_filter_id > 0;
	for (FlatArray<SeedHit>::DataIterator hit = begin; hit < end; ++hit) {
		const auto f = hit->frame;
		target.ungapped_score[f] = std::max(target.ungapped_score[f], hit->score);
		if (!diagonal_segments[f].empty() && diagonal_segments[f].back().diag() == hit->diag() && diagonal_segments[f].back().subject_end() >= hit->j)
			continue;
		const int8_t* cbs = ::Stats::CBS::hauser(config.comp_based_stats) ? query_cb[f].int8.data() : nullptr;
		const DiagonalSegment d = xdrop_ungapped(query_seq[f], cbs, target.seq, hit->i, hit->j, with_diag_filter);
		if (d.score > 0)
			diagonal_segments[f].push_back(d);
	}

	if (with_diag_filter) {
		const ApproxHsp h = Chaining::hamming_ext(diagonal_segments[0].begin(), diagonal_segments[0].end(), query_seq[0].length(), target.seq.length());
		if (h.score > 0) {
			target.done = true;
			target.hsp[0].push_back(h);
			return target;
		}
		if (h.score < 0)
			return target;
	}

	for (int frame = 0; frame < align_mode.query_contexts; ++frame) {
		if (diagonal_segments[frame].empty())
			continue;
		std::stable_sort(diagonal_segments[frame].begin(), diagonal_segments[frame].end(), DiagonalSegment::cmp_diag);
		tie(std::ignore, target.hsp[frame]) = Chaining::run(query_seq[frame], target.seq, diagonal_segments[frame].begin(), diagonal_segments[frame].end(), config.log_extend, frame);
		target.hsp[frame].sort(ApproxHsp::cmp_diag);
	}
	return target;
}

void ungapped_stage_worker(size_t i, size_t thread_id, const Sequence *query_seq, const Bias_correction *query_cb, const ::Stats::Composition* query_comp, FlatArray<SeedHit>::Iterator seed_hits, vector<uint32_t>::const_iterator target_block_ids, vector<WorkTarget> *out, mutex *mtx, Statistics* stat, const Block* targets, const Mode mode) {
	Statistics stats;
	const int16_t* query_matrix = nullptr;
	WorkTarget target = ungapped_stage(seed_hits.begin(i), seed_hits.end(i), query_seq, query_cb, *query_comp, &query_matrix, target_block_ids[i], stats, *targets, mode);
	{
		std::lock_guard<mutex> guard(*mtx);
		out->push_back(std::move(target));
		*stat += stats;
	}
	delete[] query_matrix;
}

vector<WorkTarget> ungapped_stage(const Sequence *query_seq, const Bias_correction *query_cb, const ::Stats::Composition& query_comp, FlatArray<SeedHit>::Iterator seed_hits, FlatArray<SeedHit>::Iterator seed_hits_end, vector<uint32_t>::const_iterator target_block_ids, DP::Flags flags, Statistics& stat, const Block& target_block, const Mode mode) {
	vector<WorkTarget> targets;
	const int64_t n = seed_hits_end - seed_hits;
	if(n == 0)
		return targets;
	targets.reserve(n);
	const int16_t* query_matrix = nullptr;
	if (flag_any(flags, DP::Flags::PARALLEL)) {
		mutex mtx;
		Util::Parallel::scheduled_thread_pool_auto(config.threads_, n, ungapped_stage_worker, query_seq, query_cb, &query_comp, seed_hits, target_block_ids, &targets, &mtx, &stat, &target_block, mode);
	}
	else {
		for (int64_t i = 0; i < n; ++i) {
			targets.push_back(ungapped_stage(seed_hits.begin(i), seed_hits.end(i), query_seq, query_cb, query_comp, &query_matrix, target_block_ids[i], stat, target_block, mode));
			for (const ApproxHsp& hsp : targets.back().hsp[0]) {
				Geo::assert_diag_bounds(hsp.d_max, query_seq[0].length(), targets.back().seq.length());
				Geo::assert_diag_bounds(hsp.d_min, query_seq[0].length(), targets.back().seq.length());
				assert(hsp.score > 0);
				assert(hsp.max_diag.score > 0);
			}
		}
	}

	delete[] query_matrix;
	return targets;
}

}
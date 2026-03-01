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
#include <utility>
#include <atomic>
#include <mutex>
#include "basic/config.h"
#include "stats/hauser_correction.h"
#include "target.h"
#include "dp/ungapped.h"
#include "target.h"
#include "util/parallel/thread_pool.h"
#include "chaining/chaining.h"
#include "def.h"
#include "util/geo/geo.h"

using std::array;
using std::vector;
using std::list;
using std::atomic;
using std::mutex;
using std::pair;

namespace Extension {

WorkTarget::WorkTarget(BlockId block_id, const Sequence& seq, Sequence query, Loc query_len_true_aa, const ::Stats::Composition& query_comp, Loc max_target_len, Statistics& stats, std::pmr::monotonic_buffer_resource& pool) :
	block_id(block_id),
	seq(seq),
	done(false)
{
	ungapped_score.fill(0);
	Stats::EMatrixAdjustRule rule;
	if (!config.anchored_swipe && (rule = ::Stats::adjust_matrix(query_comp, query_len_true_aa, config.comp_based_stats, seq)) != Stats::eDontAdjustMatrix) {
		matrix.reset(new ::Stats::TargetMatrix(query_comp, query_len_true_aa, config.comp_based_stats, seq, stats, pool, rule));
		/*if (config.anchored_swipe) {
			TaskTimer timer;
			profile = DP::make_profile16(query, *matrix, query.length() + max_target_len + 32);
			profile_rev = profile.reverse();
			stats.inc(Statistics::TIME_PROFILE_GENERATION, timer.microseconds());
		}*/
	}
}

WorkTarget ungapped_stage(FlatArray<SeedHit>::DataIterator begin, FlatArray<SeedHit>::DataIterator end, const Sequence *query_seq, const HauserCorrection *query_cb, const ::Stats::Composition& query_comp, uint32_t block_id, Loc max_target_len, Statistics& stats, const Block& targets, const Mode mode, std::pmr::monotonic_buffer_resource& pool, const Search::Config& cfg) {
	array<vector<DiagonalSegment>, MAX_CONTEXT> diagonal_segments;
	const SequenceSet& ref_seqs = targets.seqs(), &ref_seqs_unmasked = targets.unmasked_seqs();
	//const bool masking = config.comp_based_stats == ::Stats::CBS::COMP_BASED_STATS_AND_MATRIX_ADJUST ? ::Stats::use_seg_masking(query_seq[0], ref_seqs_unmasked[block_id]) : true;
	const bool masking = true;
	const bool with_diag_filter = (config.hamming_ext || config.diag_filter_cov.present() || config.diag_filter_id.present()) && !config.mutual_cover.present() && align_mode.query_contexts == 1;
	WorkTarget target(block_id, masking ? ref_seqs[block_id] : ref_seqs_unmasked[block_id], *query_seq, ::Stats::count_true_aa(query_seq[0]), query_comp, max_target_len, stats, pool);
	
	if (mode == Mode::FULL) {
		for (FlatArray<SeedHit>::DataIterator hit = begin; hit < end; ++hit)
			target.ungapped_score[hit->frame] = std::max(target.ungapped_score[hit->frame], hit->score);
		if (!with_diag_filter)
			return target;
	}
	if (end - begin == 1 && align_mode.query_translated) { // ???
		target.ungapped_score[begin->frame] = begin->score;
		target.hsp[begin->frame].emplace_back(begin->diag(), begin->diag(), begin->score, begin->frame, begin->query_range(), begin->target_range(), begin->diag_segment());
		return target;
	}
	std::sort(begin, end);
	const bool use_hauser = ::Stats::CBS::hauser(config.comp_based_stats);
	for (FlatArray<SeedHit>::DataIterator hit = begin; hit < end; ++hit) {
		const auto f = hit->frame;
		target.ungapped_score[f] = std::max(target.ungapped_score[f], hit->score);
		if (!diagonal_segments[f].empty() && diagonal_segments[f].back().diag() == hit->diag() && diagonal_segments[f].back().subject_end() >= hit->j)
			continue;
		const int8_t* cbs = use_hauser ? query_cb[f].int8.data() : nullptr;
		const DiagonalSegment d = xdrop_ungapped(query_seq[f], cbs, target.seq, hit->i, hit->j, with_diag_filter);
		if (d.score > 0)
			diagonal_segments[f].push_back(d);
	}

	if (with_diag_filter) {
		const ApproxHsp h = Chaining::hamming_ext(diagonal_segments[0].begin(), diagonal_segments[0].end(), query_seq[0].length(), target.seq.length(), !cfg.lin_stage1_target && !config.lin_stage1_query);
		if (h.score > 0) {
			target.done = true;
			target.hsp[0].push_back(h);
			return target;
		}
		if (h.score < 0) {
			target.ungapped_score[0] = 0;
			return target;
		}
	}

	if (mode == Mode::FULL)
		return target;

	for (int frame = 0; frame < align_mode.query_contexts; ++frame) {
		if (diagonal_segments[frame].empty())
			continue;
		std::stable_sort(diagonal_segments[frame].begin(), diagonal_segments[frame].end(), DiagonalSegment::cmp_diag);
		tie(std::ignore, target.hsp[frame]) = Chaining::run(query_seq[frame], target.seq, diagonal_segments[frame].begin(), diagonal_segments[frame].end(), config.log_extend, frame);
		target.hsp[frame].sort(ApproxHsp::cmp_diag);
	}
	return target;
}

void ungapped_stage_worker(size_t i, size_t thread_id, const Sequence *query_seq, const HauserCorrection *query_cb, const ::Stats::Composition* query_comp, FlatArray<SeedHit>::Iterator seed_hits, vector<uint32_t>::const_iterator target_block_ids, Loc max_target_len, vector<WorkTarget> *out, mutex *mtx, Statistics* stat, const Block* targets, const Mode mode, std::pmr::monotonic_buffer_resource* pool, const Search::Config *cfg) {
	Statistics stats;
	WorkTarget target = ungapped_stage(seed_hits.begin(i), seed_hits.end(i), query_seq, query_cb, *query_comp, target_block_ids[i], max_target_len, stats, *targets, mode, *pool, *cfg);
	{
		std::lock_guard<mutex> guard(*mtx);
		out->push_back(std::move(target));
		*stat += stats;
	}
}

vector<WorkTarget> ungapped_stage(const Sequence *query_seq, const HauserCorrection *query_cb, const ::Stats::Composition& query_comp, FlatArray<SeedHit>::Iterator seed_hits, FlatArray<SeedHit>::Iterator seed_hits_end, vector<uint32_t>::const_iterator target_block_ids, DP::Flags flags, Statistics& stat, const Block& target_block, const Mode mode, std::pmr::monotonic_buffer_resource& pool, const Search::Config &cfg) {
	vector<WorkTarget> targets;
	const int64_t n = seed_hits_end - seed_hits;
	if(n == 0)
		return targets;
	Loc max_target_len = 0;
	/*if (config.anchored_swipe && config.comp_based_stats == Stats::CBS::COMP_BASED_STATS_AND_MATRIX_ADJUST) {
		for(int64_t i= 0; i < n; ++i)
			max_target_len = std::max(max_target_len, target_block.seqs()[target_block_ids[i]].length());
	}*/
	targets.reserve(n);
	if (flag_any(flags, DP::Flags::PARALLEL)) {
		mutex mtx;
		Util::Parallel::scheduled_thread_pool_auto(config.threads_, n, ungapped_stage_worker, query_seq, query_cb, &query_comp, seed_hits, target_block_ids, max_target_len, &targets, &mtx, &stat, &target_block, mode, &pool, &cfg);
	}
	else {
		for (int64_t i = 0; i < n; ++i) {
			/*const double len_ratio = query_seq->length_ratio(target_block.seqs()[target_block_ids[i]]);
			if (len_ratio < config.min_length_ratio)
				continue;*/
			targets.push_back(ungapped_stage(seed_hits.begin(i), seed_hits.end(i), query_seq, query_cb, query_comp, target_block_ids[i], max_target_len, stat, target_block, mode, pool, cfg));
			for (const ApproxHsp& hsp : targets.back().hsp[0]) {
				Geo::assert_diag_bounds(hsp.d_max, query_seq[0].length(), targets.back().seq.length());
				Geo::assert_diag_bounds(hsp.d_min, query_seq[0].length(), targets.back().seq.length());
				assert(hsp.score > 0);
				assert(hsp.max_diag.score > 0);
			}
		}
	}

	return targets;
}

}
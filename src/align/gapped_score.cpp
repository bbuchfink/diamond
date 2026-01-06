/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

#include <map>
#include <numeric>
#include "target.h"
#include "dp/dp.h"
#include "util/geo/interval.h"
#include "def.h"
#include "util/geo/geo.h"
#include "stats/cbs.h"

using std::runtime_error;
using std::vector;
using std::array;
using std::list;
using std::map;
using std::endl;
using std::unique_ptr;
using std::pair;
using std::numeric_limits;

namespace Extension {

int band(int len, const Mode mode) {
	if (config.padding > 0)
		return config.padding;
	if (mode == Mode::BANDED_FAST) {
		if (len < 50)
			return 12;
		if (len < 100)
			return 16;
		if (len < 250)
			return 30;
		if (len < 350)
			return 40;
		return 64;
	}
	else {
		if (len < 50)
			return 15;
		if (len < 100)
			return 20;
		if (len < 150)
			return 30;
		if (len < 200)
			return 50;
		if (len < 250)
			return 60;
		if (len < 350)
			return 100;
		if (len < 500)
			return 120;
		return 150;
	}
}

static int hsp_band(int base_band, int qlen, int tlen, const ApproxHsp& hsp) {
	/*if (!config.classic_band) {
		return std::max(Loc((hsp.d_max - hsp.d_min) * 0.15), 32);
	}*/
	if (config.narrow_band_cov == 0.0)
		return base_band;
	if ((double)hsp.query_range.length() / qlen >= config.narrow_band_cov
		|| (double)hsp.subject_range.length() / tlen >= config.narrow_band_cov) {
		return (Loc)((hsp.d_max - hsp.d_min) * config.narrow_band_factor);
	}
	return base_band;
}

Match::Match(BlockId target_block_id, const Sequence& seq, std::unique_ptr<::Stats::TargetMatrix>&& matrix, std::array<std::list<Hsp>, MAX_CONTEXT>& hsps, int ungapped_score) :
	target_block_id(target_block_id),
	seq(seq),
	matrix(std::move(matrix)),
	filter_score(0),
	filter_evalue(DBL_MAX),
	ungapped_score(ungapped_score)
{
	if (config.max_hsps != 1)
		throw runtime_error("Match::Match max_hsps != 1.");
	for (int i = 0; i < align_mode.query_contexts; ++i)
		hsp.splice(hsp.end(), hsps[i]);
	if (hsp.empty())
		throw runtime_error("Match::Match hsp.empty()");
	hsp.sort();
	hsp.resize(1);
	filter_evalue = hsp.front().evalue;
	filter_score = hsp.front().score;
}

static void add_dp_targets(const WorkTarget& target,
	BlockId target_idx,
	const ::Stats::TargetMatrix* matrix,
	const Sequence* query_seq,
	array<DP::Targets, MAX_CONTEXT>& dp_targets,
	DP::Flags flags, const HspValues hsp_values,
	const Mode mode,
	const Search::Config& cfg)
{
	const Loc band = Extension::band(query_seq->length(), mode),
		slen = target.seq.length();
	const unsigned score_width = matrix ? matrix->score_width() : 0;
	for (int frame = 0; frame < align_mode.query_contexts; ++frame) {

		const Loc qlen = query_seq[frame].length();

		if (mode == Mode::FULL) {
			if (target.ungapped_score[frame] == 0)
				continue;
			const unsigned b = DP::BandedSwipe::bin(hsp_values, qlen, 0, target.ungapped_score[frame], (int64_t)qlen * (int64_t)slen, score_width, 0);
			//dp_targets[frame][b].emplace_back(target.seq, slen, 0, 0, Interval(), 0, target_idx, qlen, matrix);
			dp_targets[frame][b].emplace_back(target.seq, slen, 0, 0, target_idx, qlen, matrix);
			continue;
		}
		if (target.hsp[frame].empty())
			continue;
		
		int d0 = INT_MAX, d1 = INT_MIN, score = 0;
		Loc j_min = numeric_limits<Loc>::max(), j_max = numeric_limits<Loc>::min();
		Anchor anchor;

		for (const ApproxHsp &hsp : target.hsp[frame]) {
			Geo::assert_diag_bounds(hsp.d_max, qlen, slen);
			Geo::assert_diag_bounds(hsp.d_min, qlen, slen);
			assert(hsp.score > 0);
			assert(hsp.max_diag.score > 0);
			const int b = hsp_band(band, qlen, slen, hsp);
			const int b0 = std::max(hsp.d_min - b, -(slen - 1)),
				b1 = std::min(hsp.d_max + 1 + b, qlen);
			const double overlap = intersect(Interval(d0, d1), Interval(b0, b1)).length();
			if ((overlap / (d1 - d0) > config.min_band_overlap || overlap / (b1 - b0) > config.min_band_overlap) && !config.anchored_swipe) {
				d0 = std::min(d0, b0);
				d1 = std::max(d1, b1);
				score = std::max(score, hsp.score);
				j_min = std::min(j_min, hsp.subject_range.begin_);
				j_max = std::max(j_max, hsp.subject_range.end_);
				if (hsp.max_diag.score > anchor.score) anchor = hsp.max_diag;
			}
			else {
				if (d0 != INT_MAX) {
					const int64_t dp_size = (int64_t)DpTarget::banded_cols(qlen, slen, d0, d1) * int64_t(d1 - d0);
					const auto bin = DP::BandedSwipe::bin(hsp_values, d1 - d0, 0, score, dp_size, score_width, 0);
					//dp_targets[frame][bin].emplace_back(target.seq, slen, d0, d1, Interval(j_min, j_max), score, target_idx, qlen, matrix, DpTarget::CarryOver(), anchor);
					dp_targets[frame][bin].emplace_back(target.seq, slen, d0, d1, target_idx, qlen, matrix, DpTarget::CarryOver(), anchor);
					dp_targets[frame][bin].back().prof = &target.profile;
					dp_targets[frame][bin].back().prof_reverse = &target.profile_rev;
				}
				d0 = b0;
				d1 = b1;
				score = hsp.score;
				j_min = hsp.subject_range.begin_;
				j_max = hsp.subject_range.end_;
				anchor = hsp.max_diag;
			}
		}

		const int64_t dp_size = (int64_t)DpTarget::banded_cols(qlen, slen, d0, d1) * int64_t(d1 - d0);
		const auto bin = DP::BandedSwipe::bin(hsp_values, d1 - d0, 0, score, dp_size, score_width, 0);
		//dp_targets[frame][bin].emplace_back(target.seq, slen, d0, d1, Interval(j_min, j_max), score, target_idx, qlen, matrix, DpTarget::CarryOver(), anchor);
		dp_targets[frame][bin].emplace_back(target.seq, slen, d0, d1, target_idx, qlen, matrix, DpTarget::CarryOver(), anchor);
		dp_targets[frame][bin].back().prof = &target.profile;
		dp_targets[frame][bin].back().prof_reverse = &target.profile_rev;
	}
}

vector<Target> align(vector<WorkTarget>& targets, const Sequence* query_seq, const char* query_id, const HauserCorrection* query_cb, int source_query_len, DP::Flags flags, const HspValues hsp_values, const Mode mode, ThreadPool& tp, const Search::Config& cfg, Statistics& stat, std::pmr::monotonic_buffer_resource& pool) {
	array<DP::Targets, MAX_CONTEXT> dp_targets;
	vector<Target> r;
	if (targets.empty())
		return r;
	r.reserve(targets.size());
	size_t cbs_targets = 0;

	for (size_t i = 0; i < targets.size(); ++i) {
		r.emplace_back(targets[i].block_id, targets[i].seq, targets[i].ungapped_score.front(), std::move(targets[i].matrix));
		if (targets[i].done) {
			assert(targets[i].hsp[0].size() == 1);
			assert(align_mode.query_contexts == 1);
			r.back().add_hit(targets[i].hsp[0].front(), query_seq[targets[i].hsp[0].front().frame].length());
		}
		else
			add_dp_targets(targets[i], (BlockId)i, r.back().matrix.get(), query_seq, dp_targets, flags, hsp_values, mode, cfg);
		if (targets[i].matrix)
			++cbs_targets;
	}
	stat.inc(Statistics::TARGET_HITS3_CBS, cbs_targets);

	if (mode == Mode::FULL)
		flags |= DP::Flags::FULL_MATRIX;
	if (mode == Mode::GLOBAL)
		flags |= DP::Flags::SEMI_GLOBAL;

	for (int frame = 0; frame < align_mode.query_contexts; ++frame) {
		const int64_t n = std::accumulate(dp_targets[frame].begin(), dp_targets[frame].end(), (int64_t)0, [](int64_t n, const DP::TargetVec& v) { return n + v.size(); });
		if (n == 0)
			continue;
		DP::Params params{
			query_seq[frame],
			query_id,
			Frame(frame),
			source_query_len,
			::Stats::CBS::hauser(config.comp_based_stats) ? query_cb[frame].int8.data() : nullptr,
			flags,
			false,
			0,
			-1,
			hsp_values,
			stat,
			&tp
		};
		DP::AnchoredSwipe::Config acfg{ query_seq[frame],
			::Stats::CBS::hauser(config.comp_based_stats) ? query_cb[frame].int8.data() : nullptr,
			0, stat, &tp, config.comp_based_stats == Stats::CBS::COMP_BASED_STATS_AND_MATRIX_ADJUST, cfg.extension_mode, false };
		list<Hsp> hsp = config.anchored_swipe
			? DP::BandedSwipe::anchored_swipe(dp_targets[frame], acfg, pool) : DP::BandedSwipe::swipe(dp_targets[frame], params);
		while (!hsp.empty())
			r[hsp.front().swipe_target].add_hit(hsp, hsp.begin());
	}
	
	vector<Target> r2;
	r2.reserve(r.size());
	for (vector<Target>::iterator i = r.begin(); i != r.end(); ++i)
		if (i->filter_evalue != DBL_MAX) {
			if(config.max_hsps == 1 || have_coords(hsp_values))
				i->inner_culling();
			r2.push_back(std::move(*i));
		}

	return r2;
}

}
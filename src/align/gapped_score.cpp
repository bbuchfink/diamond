/****
DIAMOND protein aligner
Copyright (C) 2013-2021 Max Planck Society for the Advancement of Science e.V.
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

#include <map>
#include <numeric>
#include <algorithm>
#include "target.h"
#include "../dp/dp.h"
#include "../util/geo/interval.h"
#include "../data/reference.h"
#include "../output/output_format.h"
#include "def.h"
#include "../dp/pfscan/pfscan.h"
#include "../util/geo/geo.h"
#include "../stats/cbs.h"
#include "../dp/ungapped.h"

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
	if (config.prefix_scan && !config.classic_band) {
		return std::max(Loc((hsp.d_max - hsp.d_min) * 0.15), 32);
	}
	if (config.narrow_band_cov == 0.0)
		return base_band;
	if ((double)hsp.query_range.length() / qlen >= config.narrow_band_cov
		|| (double)hsp.subject_range.length() / tlen >= config.narrow_band_cov) {
		return (Loc)((hsp.d_max - hsp.d_min) * config.narrow_band_factor);
	}
	return base_band;
}

Match::Match(BlockId target_block_id, const Sequence& seq, const ::Stats::TargetMatrix& matrix, std::array<std::list<Hsp>, MAX_CONTEXT> &hsps, int ungapped_score):
	target_block_id(target_block_id),
	seq(seq),
	matrix(matrix),
	filter_score(0),
	filter_evalue(DBL_MAX),
	ungapped_score(ungapped_score)
{
	if (config.max_hsps != 1)
		throw std::runtime_error("Match::Match max_hsps != 1.");
	for (int i = 0; i < align_mode.query_contexts; ++i)
		hsp.splice(hsp.end(), hsps[i]);
	if (hsp.empty())
		throw std::runtime_error("Match::Match hsp.empty()");
	hsp.sort();
	hsp.resize(1);
	filter_evalue = hsp.front().evalue;
	filter_score = hsp.front().score;
}

static void add_dp_targets(const WorkTarget& target,
	int target_idx,
	const Sequence* query_seq,
	array<DP::Targets, MAX_CONTEXT>& dp_targets,
	DP::Flags flags, const HspValues hsp_values,
	const Mode mode,
	const Search::Config& cfg)
{
	const Loc band = Extension::band(query_seq->length(), mode),
		slen = target.seq.length();
	const ::Stats::TargetMatrix* matrix = target.matrix.scores.empty() ? nullptr : &target.matrix;
	const unsigned score_width = matrix ? matrix->score_width() : 0;
	for (int frame = 0; frame < align_mode.query_contexts; ++frame) {

		const Loc qlen = query_seq[frame].length();

		if (mode == Mode::FULL) {
			if (target.ungapped_score[frame] == 0)
				continue;
			const unsigned b = DP::BandedSwipe::bin(hsp_values, qlen, 0, target.ungapped_score[frame], (int64_t)qlen * (int64_t)slen, score_width, 0);
			dp_targets[frame][b].emplace_back(target.seq, slen, 0, 0, Interval(), 0, target_idx, qlen, matrix);
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
					dp_targets[frame][bin].emplace_back(target.seq, slen, d0, d1, Interval(j_min, j_max), score, target_idx, qlen, matrix, DpTarget::CarryOver(), anchor);
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
		dp_targets[frame][bin].emplace_back(target.seq, slen, d0, d1, Interval(j_min, j_max), score, target_idx, qlen, matrix, DpTarget::CarryOver(), anchor);
	}
}

pair<vector<Target>, Stats> align(const vector<WorkTarget> &targets, const Sequence *query_seq, const char* query_id, const Bias_correction *query_cb, int source_query_len, DP::Flags flags, const HspValues hsp_values, const Mode mode, ThreadPool& tp, const Search::Config& cfg, Statistics &stat) {
	array<DP::Targets, MAX_CONTEXT> dp_targets;
	vector<Target> r;
	if (targets.empty())
		return make_pair(r, Stats());
	r.reserve(targets.size());
	size_t cbs_targets = 0;

	for (int i = 0; i < (int)targets.size(); ++i) {
		r.emplace_back(targets[i].block_id, targets[i].seq, targets[i].ungapped_score.front(), targets[i].matrix);
		if (targets[i].done)
			r.back().add_hit(targets[i].hsp[0].front(), query_seq[targets[i].hsp[0].front().frame].length());
		else
			add_dp_targets(targets[i], i, query_seq, dp_targets, flags, hsp_values, mode, cfg);
		if (targets[i].adjusted_matrix())
			++cbs_targets;		
	}
	stat.inc(Statistics::TARGET_HITS3_CBS, cbs_targets);

	if (mode == Mode::FULL)
		flags |= DP::Flags::FULL_MATRIX;
	if (mode == Mode::GLOBAL)
		flags |= DP::Flags::SEMI_GLOBAL;

	Stats stats;

	for (int frame = 0; frame < align_mode.query_contexts; ++frame) {
		const int64_t n = std::accumulate(dp_targets[frame].begin(), dp_targets[frame].end(), (int64_t)0, [](int64_t n, const vector<DpTarget>& v) { return n + v.size(); });
		if (n == 0)
			continue;
		stats.extension_count += n;
		if (config.prefix_scan) {
			LongScoreProfile<int16_t> p;
			LongScoreProfile<int8_t> p8;
			const bool hauser_cbs = ::Stats::CBS::hauser(config.comp_based_stats);
			TaskTimer timer;
			p = DP::make_profile16(query_seq[frame], hauser_cbs ? query_cb[frame].int8.data() : nullptr, 0);
			p8 = DP::make_profile8(query_seq[frame], hauser_cbs ? query_cb[frame].int8.data() : nullptr, 0);
			LongScoreProfile<int16_t> p_rev(p.reverse());
			LongScoreProfile<int8_t> p8_rev(p8.reverse());
			stat.inc(Statistics::TIME_PROFILE, timer.microseconds());
			const auto v = p.pointers(0), vr = p_rev.pointers(0);
			const auto v8 = p8.pointers(0), vr8 = p8_rev.pointers(0);
			for (int i = 0; i < 6; ++i)
				for (const DpTarget& t : dp_targets[frame][i]) {
					const char* tid = cfg.target->ids()[r[t.target_idx].block_id];
					DP::PrefixScan::Config cfg{ query_seq[frame], t.seq, query_id, tid, t.d_begin, t.d_end,
						t.chaining_target_range, v.data(), vr.data(), v8.data(), vr8.data(), stat, 0, 0, t.chaining_score };
					//Hsp h = DP::PrefixScan::align(cfg);
					//const DiagonalSegment anchor = make_null_anchor(t.anchor);
					const Anchor anchor = make_clipped_anchor(t.anchor, query_seq[frame], hauser_cbs ? query_cb[frame].int8.data() : nullptr, t.seq);
					if (anchor.score == 0)
						continue;
					Hsp h = DP::PrefixScan::align_anchored(anchor, cfg);
					if (h.evalue <= config.max_evalue)
						r[t.target_idx].add_hit(std::move(h));
				}
		}
		else {
			DP::Params params{
				query_seq[frame],
				query_id,
				Frame(frame),
				source_query_len,
				::Stats::CBS::hauser(config.comp_based_stats) ? query_cb[frame].int8.data() : nullptr,
				flags,
				hsp_values,
				stat,
				&tp
			};
			DP::AnchoredSwipe::Config cfg{ query_seq[frame], ::Stats::CBS::hauser(config.comp_based_stats) ? query_cb[frame].int8.data() : nullptr, 0, stat, &tp };
			list<Hsp> hsp = config.anchored_swipe ? DP::BandedSwipe::anchored_swipe(dp_targets[frame], cfg) : DP::BandedSwipe::swipe(dp_targets[frame], params);
			while (!hsp.empty())
				r[hsp.front().swipe_target].add_hit(hsp, hsp.begin());
		}
	}

	vector<Target> r2;
	r2.reserve(r.size());
	for (vector<Target>::iterator i = r.begin(); i != r.end(); ++i)
		if (i->filter_evalue != DBL_MAX) {
			if(config.max_hsps == 1 || have_coords(hsp_values))
				i->inner_culling();
			r2.push_back(std::move(*i));
		}

	return make_pair(r2, stats);
}

}
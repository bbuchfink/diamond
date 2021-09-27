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
#include <algorithm>
#include "target.h"
#include "../dp/dp.h"
#include "../util/interval.h"
#include "../data/reference.h"
#include "../output/output_format.h"

using std::vector;
using std::array;
using std::list;
using std::map;
using std::endl;

namespace Extension {

static int band(int len, const Mode mode) {
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

Match::Match(size_t target_block_id, const Sequence& seq, const Stats::TargetMatrix& matrix, std::array<std::list<Hsp>, MAX_CONTEXT> &hsps, int ungapped_score):
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

static void add_dp_targets(const WorkTarget& target, int target_idx, const Sequence* query_seq, array<DP::Targets, MAX_CONTEXT>& dp_targets, DP::Flags flags, const HspValues hsp_values, const Mode mode) {
	const int band = Extension::band((int)query_seq->length(), mode),
		slen = (int)target.seq.length();
	const size_t qlen = query_seq[0].length();
	const Stats::TargetMatrix* matrix = target.matrix.scores.empty() ? nullptr : &target.matrix;
	const unsigned score_width = matrix ? matrix->score_width() : 0;
	for (int frame = 0; frame < align_mode.query_contexts; ++frame) {

		if (mode == Mode::FULL) {
			if (target.ungapped_score[frame] == 0)
				continue;
			const unsigned b = DP::BandedSwipe::bin(hsp_values, query_seq[0].length(), 0, target.ungapped_score[frame], query_seq[0].length() * target.seq.length(), score_width, 0);
			dp_targets[frame][b].emplace_back(target.seq, target.seq.length(), 0, 0, target_idx, (int)query_seq->length(), matrix);
			continue;
		}
		if (target.hsp[frame].empty())
			continue;
		const int qlen = (int)query_seq[frame].length();
		int d0 = INT_MAX, d1 = INT_MIN, j0 = INT_MAX, j1 = INT_MIN, bits = 0;

		for (const Hsp_traits &hsp : target.hsp[frame]) {
			const int b0 = std::max(hsp.d_min - band, -(slen - 1)),
				b1 = std::min(hsp.d_max + 1 + band, qlen);
			const double overlap = intersect(interval(d0, d1), interval(b0, b1)).length();
			if (overlap / (d1 - d0) > config.min_band_overlap || overlap / (b1 - b0) > config.min_band_overlap) {
				d0 = std::min(d0, b0);
				d1 = std::max(d1, b1);
				j0 = std::min(j0, hsp.subject_range.begin_);
				j1 = std::max(j1, hsp.subject_range.end_);
				bits = std::max(bits, (int)DP::BandedSwipe::bin(hsp_values, d1 - d0, 0, hsp.score, size_t(d1 - d0) * target.seq.length(), score_width, 0));
			}
			else {
				if (d0 != INT_MAX)
					dp_targets[frame][bits].emplace_back(target.seq, target.seq.length(), d0, d1, target_idx, (int)query_seq->length(), matrix);
				d0 = b0;
				d1 = b1;
				j0 = hsp.subject_range.begin_;
				j1 = hsp.subject_range.end_;
				bits = (int)DP::BandedSwipe::bin(hsp_values, d1 - d0, 0, hsp.score, size_t(d1 - d0) * target.seq.length(), score_width, 0);
			}
		}

		dp_targets[frame][bits].emplace_back(target.seq, target.seq.length(), d0, d1, target_idx, (int)query_seq->length(), matrix);
	}
}

vector<Target> align(const vector<WorkTarget> &targets, const Sequence *query_seq, const Bias_correction *query_cb, int source_query_len, DP::Flags flags, const HspValues hsp_values, const Mode mode, Statistics &stat) {
	array<DP::Targets, MAX_CONTEXT> dp_targets;
	vector<Target> r;
	if (targets.empty())
		return r;
	r.reserve(targets.size());
	size_t cbs_targets = 0;
	for (int i = 0; i < (int)targets.size(); ++i) {
		add_dp_targets(targets[i], i, query_seq, dp_targets, flags, hsp_values, mode);
		if (targets[i].adjusted_matrix())
			++cbs_targets;
		r.emplace_back(targets[i].block_id, targets[i].seq, targets[i].ungapped_score.front(), targets[i].matrix);
	}
	stat.inc(Statistics::TARGET_HITS3_CBS, cbs_targets);

	if (mode == Mode::FULL)
		flags |= DP::Flags::FULL_MATRIX;
	if (mode == Mode::GLOBAL)
		flags |= DP::Flags::SEMI_GLOBAL;

	for (int frame = 0; frame < align_mode.query_contexts; ++frame) {
		if (dp_targets[frame].empty())
			continue;
		DP::Params params{
			query_seq[frame],
			Frame(frame),
			source_query_len,
			Stats::CBS::hauser(config.comp_based_stats) ? query_cb[frame].int8.data() : nullptr,
			flags,
			hsp_values,
			stat
		};
		list<Hsp> hsp = DP::BandedSwipe::swipe(dp_targets[frame], params);
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
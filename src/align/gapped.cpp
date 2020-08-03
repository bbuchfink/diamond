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

#include <algorithm>
#include "target.h"
#include "../dp/dp.h"
#include "../util/interval.h"
#include "../data/reference.h"

using std::vector;
using std::array;
using std::list;

namespace Extension {

int band(int len) {
	if (config.padding > 0)
		return config.padding;
	if ((config.sensitivity <= Sensitivity::SENSITIVE && config.ext != "banded-slow") || config.ext == "banded-fast") {
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

Match::Match(size_t target_block_id, bool outranked, std::array<std::list<Hsp>, MAX_CONTEXT> &hsps, int ungapped_score):
	target_block_id(target_block_id),
	filter_score(0),
	ungapped_score(ungapped_score),
	outranked(outranked)
{
	for (unsigned i = 0; i < align_mode.query_contexts; ++i)
		hsp.splice(hsp.end(), hsps[i]);
	hsp.sort();
	if (!hsp.empty())
		filter_score = hsp.front().score;
	if (config.max_hsps > 0)
		max_hsp_culling();
}

void add_dp_targets(const WorkTarget &target, int target_idx, const sequence *query_seq, array<array<vector<DpTarget>, 2>, MAX_CONTEXT> &dp_targets) {
	const int band = Extension::band((int)query_seq->length()),
		slen = (int)target.seq.length();
	for (unsigned frame = 0; frame < align_mode.query_contexts; ++frame) {

		if (config.ext == "full" && (int)query_seq[0].length() < config.full_sw_len) {
			dp_targets[frame][0].emplace_back(target.seq, 0, 0, 0, 0, target_idx, (int)query_seq->length());
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
			if (overlap / (d1 - d0) >= config.min_band_overlap || overlap / (b1 - b0) >= config.min_band_overlap) {
				d0 = std::min(d0, b0);
				d1 = std::max(d1, b1);
				j0 = std::min(j0, hsp.subject_range.begin_);
				j1 = std::max(j1, hsp.subject_range.end_);
				bits = std::max(bits, hsp.score <= config.cutoff_score_8bit ? 0 : 1);
			}
			else {
				if (d0 != INT_MAX)
					dp_targets[frame][bits].emplace_back(target.seq, d0, d1, j0, j1, target_idx, (int)query_seq->length());
				d0 = b0;
				d1 = b1;
				j0 = hsp.subject_range.begin_;
				j1 = hsp.subject_range.end_;
				bits = hsp.score <= config.cutoff_score_8bit ? 0 : 1;
			}
		}

		dp_targets[frame][bits].emplace_back(target.seq, d0, d1, j0, j1, target_idx, (int)query_seq->length());

	}
}

vector<Target> align(const vector<WorkTarget> &targets, const sequence *query_seq, const Bias_correction *query_cb, int flags, Statistics &stat) {
	const int raw_score_cutoff = score_matrix.rawscore(config.min_bit_score == 0 ? score_matrix.bitscore(config.max_evalue, (unsigned)query_seq[0].length()) : config.min_bit_score);

	array<array<vector<DpTarget>, 2>, MAX_CONTEXT> dp_targets;
	vector<Target> r;
	if (targets.empty())
		return r;
	r.reserve(targets.size());
	for (int i = 0; i < (int)targets.size(); ++i) {
		add_dp_targets(targets[i], i, query_seq, dp_targets);
		r.emplace_back(targets[i].block_id, targets[i].seq, targets[i].outranked, targets[i].ungapped_score);
	}

	if (config.ext == "full" && (int)query_seq[0].length() < config.full_sw_len)
		flags |= DP::FULL_MATRIX;

	for (unsigned frame = 0; frame < align_mode.query_contexts; ++frame) {
		if (dp_targets[frame].empty())
			continue;
		list<Hsp> hsp = DP::BandedSwipe::swipe(
			query_seq[frame],
			dp_targets[frame][0],
			dp_targets[frame][1],
			Frame(frame),
			config.comp_based_stats ? &query_cb[frame] : nullptr,
			flags,
			raw_score_cutoff,
			stat);
		while (!hsp.empty())
			r[hsp.front().swipe_target].add_hit(hsp, hsp.begin());
	}

	vector<Target> r2;
	r2.reserve(r.size());
	for (vector<Target>::iterator i = r.begin(); i != r.end(); ++i)
		if (i->filter_score > 0)
			r2.push_back(std::move(*i));

	return r2;
}

void add_dp_targets(const Target &target, int target_idx, const sequence *query_seq, array<array<vector<DpTarget>, 2>, MAX_CONTEXT> &dp_targets) {
	const int band = Extension::band((int)query_seq->length()),
		slen = (int)target.seq.length();
	for (unsigned frame = 0; frame < align_mode.query_contexts; ++frame) {
		const int qlen = (int)query_seq[frame].length();
		for (const Hsp &hsp : target.hsp[frame]) {
			vector<DpTarget>& v = hsp.score < 255 && (hsp.query_range.end_ - hsp.query_range.begin_ < 256) ? dp_targets[frame][0] : dp_targets[frame][1];
			v.emplace_back(target.seq, hsp.query_range.begin_, hsp.query_range.end_, hsp.subject_range.begin_, hsp.subject_range.end_, target_idx);
		}
	}
}

vector<Match> align(vector<Target> &targets, const sequence *query_seq, const Bias_correction *query_cb, int source_query_len, int flags, Statistics &stat) {
	const int raw_score_cutoff = score_matrix.rawscore(config.min_bit_score == 0 ? score_matrix.bitscore(config.max_evalue, (unsigned)query_seq[0].length()) : config.min_bit_score);

	array<array<vector<DpTarget>, 2>, MAX_CONTEXT> dp_targets;
	vector<Match> r;
	if (targets.empty())
		return r;
	r.reserve(targets.size());

	if (config.traceback_mode == TracebackMode::SCORE_ONLY || config.ext == "full") {
		for (Target &t : targets)
			r.emplace_back(t.block_id, t.outranked, t.hsp, t.ungapped_score);
		return r;
	}

	for (int i = 0; i < (int)targets.size(); ++i) {
		if (config.log_subject)
			std::cout << "Target=" << ref_ids::get()[targets[i].block_id] << " id=" << i << endl;
		add_dp_targets(targets[i], i, query_seq, dp_targets);
		r.emplace_back(targets[i].block_id, targets[i].outranked, targets[i].ungapped_score);
	}

	for (unsigned frame = 0; frame < align_mode.query_contexts; ++frame) {
		if (dp_targets[frame].empty())
			continue;
		list<Hsp> hsp = DP::BandedSwipe::swipe(
			query_seq[frame],
			dp_targets[frame][0],
			dp_targets[frame][1],
			Frame(frame),
			config.comp_based_stats ? &query_cb[frame] : nullptr,
			DP::TRACEBACK | flags,
			raw_score_cutoff,
			stat);
		while (!hsp.empty())
			r[hsp.front().swipe_target].add_hit(hsp, hsp.begin());
	}

	for (Match &match : r)
		match.inner_culling(source_query_len);

	return r;
}

}

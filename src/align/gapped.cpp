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

Match::Match(size_t target_block_id, std::array<std::list<Hsp>, MAX_CONTEXT> &hsps, int ungapped_score):
	target_block_id(target_block_id),
	filter_score(0),
	filter_evalue(DBL_MAX),
	ungapped_score(ungapped_score)
{
	for (unsigned i = 0; i < align_mode.query_contexts; ++i)
		hsp.splice(hsp.end(), hsps[i]);
	hsp.sort();
	if (!hsp.empty()) {
		filter_evalue = hsp.front().evalue;
		filter_score = hsp.front().score;
	}
	if (config.max_hsps > 0)
		max_hsp_culling();
}

void add_dp_targets(const WorkTarget& target, int target_idx, const Sequence* query_seq, array<array<vector<DpTarget>, 3>, MAX_CONTEXT>& dp_targets, int flags) {
	const int band = Extension::band((int)query_seq->length()),
		slen = (int)target.seq.length();
	const size_t qlen = query_seq[0].length();
	const Stats::TargetMatrix* matrix = target.matrix.scores.empty() ? nullptr : &target.matrix;
	for (unsigned frame = 0; frame < align_mode.query_contexts; ++frame) {

		if (config.ext == "full") {
			if (target.ungapped_score[frame] == 0)
				continue;
			int b = target.ungapped_score[frame] <= config.cutoff_score_8bit ? 0 : 1;
			if (((flags & DP::TRACEBACK) || (flags & DP::WITH_COORDINATES)) && query_seq[0].length() >= 256)
				b = 1;
			if ((flags & DP::TRACEBACK) && query_seq[0].length() * target.seq.length() > config.max_swipe_dp)
				b = 2;
			dp_targets[frame][b].emplace_back(target.seq, 0, 0, 0, 0, target_idx, (int)query_seq->length(), matrix);
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
				bits = std::max(bits, (hsp.score <= config.cutoff_score_8bit && (d1 - d0) < 256) ? 0 : 1);
				if (flags & DP::TRACEBACK && size_t(d1 - d0) * qlen > config.max_swipe_dp)
					bits = 2;
			}
			else {
				if (d0 != INT_MAX)
					dp_targets[frame][bits].emplace_back(target.seq, d0, d1, j0, j1, target_idx, (int)query_seq->length(), matrix);
				d0 = b0;
				d1 = b1;
				j0 = hsp.subject_range.begin_;
				j1 = hsp.subject_range.end_;
				bits = (hsp.score <= config.cutoff_score_8bit && (d1 - d0) < 256) ? 0 : 1;
				if (flags & DP::TRACEBACK && size_t(d1 - d0) * qlen > config.max_swipe_dp)
					bits = 2;
			}
		}

		dp_targets[frame][bits].emplace_back(target.seq, d0, d1, j0, j1, target_idx, (int)query_seq->length(), matrix);
	}
}

vector<Target> align(const vector<WorkTarget> &targets, const Sequence *query_seq, const Bias_correction *query_cb, int source_query_len, int flags, Statistics &stat) {
	array<array<vector<DpTarget>, 3>, MAX_CONTEXT> dp_targets;
	vector<Target> r;
	if (targets.empty())
		return r;
	r.reserve(targets.size());
	size_t cbs_targets = 0;
	for (int i = 0; i < (int)targets.size(); ++i) {
		add_dp_targets(targets[i], i, query_seq, dp_targets, flags);
		if (targets[i].adjusted_matrix())
			++cbs_targets;
		r.emplace_back(targets[i].block_id, targets[i].seq, targets[i].ungapped_score.front(), targets[i].matrix);
	}
	stat.inc(Statistics::TARGET_HITS3_CBS, cbs_targets);

	if (config.ext == "full")
		flags |= DP::FULL_MATRIX;

	for (unsigned frame = 0; frame < align_mode.query_contexts; ++frame) {
		if (dp_targets[frame].empty())
			continue;
		list<Hsp> hsp = DP::BandedSwipe::swipe(
			query_seq[frame],
			dp_targets[frame][0],
			dp_targets[frame][1],
			dp_targets[frame][2],
			nullptr,
			Frame(frame),
			Stats::CBS::hauser(config.comp_based_stats) ? &query_cb[frame] : nullptr,
			flags,
			stat);
		while (!hsp.empty())
			r[hsp.front().swipe_target].add_hit(hsp, hsp.begin());
	}

	vector<Target> r2;
	r2.reserve(r.size());
	for (vector<Target>::iterator i = r.begin(); i != r.end(); ++i)
		if (i->filter_evalue != DBL_MAX) {
			if((flags & DP::TRACEBACK) || (flags & DP::WITH_COORDINATES))
				i->inner_culling(source_query_len);
			r2.push_back(std::move(*i));
		}

	return r2;
}

vector<Target> full_db_align(const Sequence *query_seq, const Bias_correction *query_cb, int flags, Statistics &stat) {	
	vector<DpTarget> v;
	vector<Target> r;
	Stats::TargetMatrix matrix;
	list<Hsp> hsp;
	
	for (unsigned frame = 0; frame < align_mode.query_contexts; ++frame) {
		ContainerIterator<DpTarget, SequenceSet> target_it(*ref_seqs::data_, ref_seqs::data_->get_length());
		list<Hsp> frame_hsp = DP::BandedSwipe::swipe(
			query_seq[frame],
			v,
			v,
			v,
			&target_it,
			Frame(frame),
			Stats::CBS::hauser(config.comp_based_stats) ? &query_cb[frame] : nullptr,
			flags | DP::FULL_MATRIX,
			stat);
		hsp.splice(hsp.begin(), frame_hsp, frame_hsp.begin(), frame_hsp.end());
	}

	map<unsigned, unsigned> subject_idx;
	while (!hsp.empty()) {
		size_t block_id = hsp.begin()->swipe_target;
		const auto it = subject_idx.emplace(block_id, (unsigned)r.size());
		if (it.second)
			r.emplace_back(block_id, ref_seqs::get()[block_id], 0, matrix);
		unsigned i = it.first->second;
		r[i].add_hit(hsp, hsp.begin());
	}

	return r;
}

void add_dp_targets(const Target &target, int target_idx, const Sequence *query_seq, array<array<vector<DpTarget>, 3>, MAX_CONTEXT> &dp_targets, int flags) {
	const int band = Extension::band((int)query_seq->length()),
		slen = (int)target.seq.length();
	const Stats::TargetMatrix* matrix = target.adjusted_matrix() ? &target.matrix : nullptr;
	for (unsigned frame = 0; frame < align_mode.query_contexts; ++frame) {
		const int qlen = (int)query_seq[frame].length();
		for (const Hsp &hsp : target.hsp[frame]) {
			const bool byte_row_counter = config.ext == "full" ? qlen < 256 : (hsp.d_end - hsp.d_begin < 256);
			int b = (hsp.score < 255 && byte_row_counter) ? 0 : 1;
			const size_t dp_size = config.ext == "full" ? query_seq[0].length() * target.seq.length() : query_seq[0].length() * size_t(hsp.d_end - hsp.d_begin);
			if (dp_size > config.max_swipe_dp && (flags & DP::TRACEBACK))
				b = 2;
			dp_targets[frame][b].emplace_back(target.seq, hsp.d_begin, hsp.d_end, hsp.seed_hit_range.begin_, hsp.seed_hit_range.end_, target_idx, qlen, matrix);
		}
	}
}

vector<Match> align(vector<Target> &targets, const Sequence *query_seq, const Bias_correction *query_cb, int source_query_len, int flags, Statistics &stat) {
	array<array<vector<DpTarget>, 3>, MAX_CONTEXT> dp_targets;
	vector<Match> r;
	if (targets.empty())
		return r;
	r.reserve(targets.size());

	if ((output_format->hsp_values == Output::NONE && config.max_hsps == 1) || (flags & DP::TRACEBACK) || ((flags & DP::WITH_COORDINATES) && !(output_format->hsp_values & Output::STATS_OR_TRANSCRIPT))) {
		for (Target &t : targets)
			r.emplace_back(t.block_id, t.hsp, t.ungapped_score);
		return r;
	}

	if (config.ext == "full") {
		flags |= DP::FULL_MATRIX;
		if ((output_format->hsp_values & Output::TRANSCRIPT) || (output_format->hsp_values & Output::STATS))
			flags |= DP::TRACEBACK;
		else
			flags |= DP::WITH_COORDINATES;
	}
	else
		flags |= DP::TRACEBACK;

	for (int i = 0; i < (int)targets.size(); ++i) {
		if (config.log_subject)
			std::cout << "Target=" << ref_ids::get()[targets[i].block_id] << " id=" << i << endl;
		add_dp_targets(targets[i], i, query_seq, dp_targets, flags);
		r.emplace_back(targets[i].block_id, targets[i].ungapped_score);
	}

	for (unsigned frame = 0; frame < align_mode.query_contexts; ++frame) {
		if (dp_targets[frame].empty())
			continue;
		list<Hsp> hsp = DP::BandedSwipe::swipe(
			query_seq[frame],
			dp_targets[frame][0],
			dp_targets[frame][1],
			dp_targets[frame][2],
			nullptr,
			Frame(frame),
			Stats::CBS::hauser(config.comp_based_stats) ? &query_cb[frame] : nullptr,
			flags,
			stat);
		while (!hsp.empty())
			r[hsp.front().swipe_target].add_hit(hsp, hsp.begin());
	}

	for (Match &match : r)
		match.inner_culling(source_query_len);

	return r;
}

}

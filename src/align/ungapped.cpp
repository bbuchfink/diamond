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
#include "../basic/config.h"
#include "../dp/comp_based_stats.h"
#include "target.h"
#include "../dp/ungapped.h"
#include "../dp/dp.h"
#include "target.h"
#include "../data/reference.h"

using std::array;
using std::vector;
using std::list;

namespace Extension {

WorkTarget ungapped_stage(Trace_pt_list::iterator begin, Trace_pt_list::iterator end, uint64_t target_offset, const sequence *query_seq, const Bias_correction *query_cb, const sequence &target_seq) {
	array<vector<Diagonal_segment>, MAX_CONTEXT> diagonal_segments;
	WorkTarget target(target_seq);
	for (Trace_pt_list::const_iterator i = begin; i < end; ++i) {
		const unsigned frame = i->frame();
		const Diagonal_segment d = xdrop_ungapped(query_seq[frame], target_seq, i->seed_offset_, uint64_t(i->subject_) - target_offset);
		if (d.score >= config.min_ungapped_raw_score)
			diagonal_segments[frame].push_back(d);
	}
	for (int frame = 0; frame < align_mode.query_contexts; ++frame) {
		if (diagonal_segments[frame].empty())
			continue;
		std::stable_sort(diagonal_segments[frame].begin(), diagonal_segments[frame].end(), Diagonal_segment::cmp_diag);
		pair<int, list<Hsp_traits>> hsp = greedy_align(query_seq[frame], query_cb[frame], target_seq, diagonal_segments[frame].begin(), diagonal_segments[frame].end(), false, frame);
		target.filter_score = std::max(target.filter_score, hsp.first);
		target.hsp[frame] = std::move(hsp.second);
	}
	return target;
}

vector<WorkTarget> ungapped_stage(const sequence *query_seq, const Bias_correction *query_cb, Trace_pt_list::iterator begin, Trace_pt_list::iterator end) {
	vector<WorkTarget> targets;
	if (begin >= end)
		return targets;
	std::sort(begin, end, hit::cmp_subject);
	Trace_pt_list::iterator i = begin, target_begin = begin;
	size_t target = ref_seqs::data_->local_position(i->subject_).first;
	while (++i < end) {
		std::pair<size_t, size_t> l = ref_seqs::data_->local_position(i->subject_);
		if (l.first != target) {
			targets.push_back(ungapped_stage(target_begin, i, ref_seqs::data_->position(target, 0), query_seq, query_cb, ref_seqs::get()[target]));
			target = l.first;
			target_begin = i;
		}
	}
	targets.push_back(ungapped_stage(target_begin, i, ref_seqs::data_->position(target, 0), query_seq, query_cb, ref_seqs::get()[target]));
	return targets;
}

}
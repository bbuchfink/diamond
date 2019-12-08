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
#include <vector>
#include "../basic/config.h"
#include "target.h"
#include "../basic/diagonal_segment.h"
#include "../dp/ungapped.h"

using std::vector;

namespace Extension {

void Target::add_hit(unsigned context, const Diagonal_segment &d) {
	filter_score_ = std::max(filter_score_, d.score);
}

void Target::ungapped_stage(unsigned context, std::vector<Diagonal_segment>::const_iterator begin, std::vector<Diagonal_segment>::const_iterator end) {

}

Target::Target(Trace_pt_list::const_iterator begin, Trace_pt_list::const_iterator end, uint64_t target_offset, const sequence *query_seq, const sequence &target_seq, const std::set<unsigned> &taxon_rank_ids) :
	filter_score_(0),
	taxon_rank_ids(taxon_rank_ids)
{
	vector<Diagonal_segment> diagonal_segments[6];
	for (Trace_pt_list::const_iterator i = begin; i < end; ++i) {
		const unsigned frame = i->frame();
		const Diagonal_segment d = xdrop_ungapped(query_seq[frame], target_seq, i->seed_offset_, uint64_t(i->subject_) - target_offset);
		if (d.score >= config.min_ungapped_raw_score) {
			add_hit(frame, d);
			diagonal_segments[frame].push_back(d);
		}
	}
	for (unsigned i = 0; i < align_mode.query_contexts; ++i)
		if (!diagonal_segments[i].empty())
			ungapped_stage(i, diagonal_segments[i].cbegin(), diagonal_segments[i].cend());
}

}

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

#ifndef TARGET_H_
#define TARGET_H_

#include <array>
#include <set>
#include <vector>
#include <stdint.h>
#include "../search/trace_pt_buffer.h"
#include "../basic/diagonal_segment.h"

namespace Extension {

struct Context {
};

struct Target {

	Target(Trace_pt_list::const_iterator begin, Trace_pt_list::const_iterator end, uint64_t target_offset, const sequence *query_seq, const sequence &target_seq, const std::set<unsigned> &taxon_rank_ids);
	void add_hit(unsigned context, const Diagonal_segment &d);
	void ungapped_stage(unsigned context, std::vector<Diagonal_segment>::const_iterator begin, std::vector<Diagonal_segment>::const_iterator end);

private:
	int filter_score_;
	std::array<Context, 6> context_;
	const std::set<unsigned> taxon_rank_ids;

};

}

#endif
/****
DIAMOND protein aligner
Copyright (C) 2016-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink

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

#pragma once
#include <stddef.h>
#include "../util/simd.h"
#include "../dp/ungapped.h"
#include "../basic/shape_config.h"
#include "../basic/statistics.h"
#include "trace_pt_buffer.h"
#include "../util/algo/pattern_matcher.h"
#include "../util/data_structures/flat_array.h"
#include "../util/scores/cutoff_table.h"

namespace Search {

struct Context {
	const PatternMatcher previous_matcher, current_matcher;
	const Util::Scores::CutoffTable cutoff_table;
	const int short_query_ungapped_cutoff;
};

}

struct Stage1_hit
{
	Stage1_hit(unsigned q_ref, unsigned q_offset, unsigned s_ref, unsigned s_offset) :
		q(q_ref + q_offset),
		s(s_ref + s_offset)
	{}
	bool operator<(const Stage1_hit &rhs) const
	{
		return q < rhs.q;
	}
	struct Query
	{
		unsigned operator()(const Stage1_hit &x) const
		{
			return x.q;
		}
	};
	unsigned q, s;
};

void search_shape(unsigned sid, unsigned query_block, char *query_buffer, char *ref_buffer);
bool use_single_indexed(double coverage, size_t query_letters, size_t ref_letters);
void setup_search();
void setup_search_cont();

namespace Search {

DECL_DISPATCH(void, stage1, (const Packed_loc* q, size_t nq, const Packed_loc* s, size_t ns, Statistics& stats, Trace_pt_buffer::Iterator& out, const unsigned sid, const Context& context))

}

extern const double SINGLE_INDEXED_SEED_SPACE_MAX_COVERAGE;
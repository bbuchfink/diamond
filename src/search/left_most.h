/****
DIAMOND protein aligner
Copyright (C) 2020 Max Planck Society for the Advancement of Science e.V.

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

#include <algorithm>
#include "search.h"
#include "sse_dist.h"
#include "../util/sequence/sequence.h"
#include "../basic/config.h"
#include "collision.h"

namespace Search {

static inline bool left_most_filter(const sequence &query,
	const Letter* subject,
	const int seed_offset,
	const int seed_len,
	const Context &context,
	bool first_shape,
	size_t shape_id)
{
	constexpr int WINDOW_LEFT = 16, WINDOW_RIGHT = 32;

	int d = std::max(seed_offset - WINDOW_LEFT, 0), window_left = std::min(WINDOW_LEFT, seed_offset);
	const Letter *q = query.data() + d, *s = subject + d;
	int window = (int)query.length() - d;
	window = std::min(window, window_left + 1 + WINDOW_RIGHT);

	const sequence subject_clipped = Util::Sequence::clip(s, window, window_left);
	window -= s + window - subject_clipped.end();

	d = subject_clipped.data() - s;
	q += d;
	s += d;
	window_left -= d;
	window -= d;

	if (!config.beta)
		return is_primary_hit(q, s, window_left, (unsigned)shape_id, window);

	const uint64_t match_mask = reduced_match(q, s, window),
		query_seed_mask = ~seed_mask(q, window);

	const uint32_t len_left = window_left + seed_len - 1,
		match_mask_left = ((1llu << len_left) - 1) & match_mask,
		query_mask_left = ((1llu << len_left) - 1) & query_seed_mask;

	const uint32_t left_hit = context.current_matcher.hit(match_mask_left, len_left) & query_mask_left;

	if (first_shape)
		return left_hit == 0;

	const uint32_t len_right = window - window_left - 1,
		match_mask_right = match_mask >> (window_left + 1),
		query_mask_right = query_seed_mask >> (window_left + 1);
	
	const uint32_t right_hit = context.previous_matcher.hit(match_mask_right, len_right) & query_mask_right;
	return left_hit == 0 && right_hit == 0;
}

}
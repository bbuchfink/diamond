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
#include "finger_print.h"

namespace Search {

static inline bool verify_hit(const Letter* q, const Letter* s, int score_cutoff, bool left, uint32_t match_mask, int sid, bool chunked, unsigned hamming_filter_id) {
	if (chunked) {
		if ((shapes[sid].mask_ & match_mask) == shapes[sid].mask_) {
			PackedSeed seed;
			if (!shapes[sid].set_seed(seed, s))
				return false;
			if (left && !current_range.lower_or_equal(seed_partition(seed)))
				return false;
			if (!left && !current_range.lower(seed_partition(seed)))
				return false;
		}
	}
	FingerPrint fq(q), fs(s);
	const unsigned id = fq.match(fs);
	return id >= hamming_filter_id;
}

static inline bool verify_hits(uint32_t mask, const Letter* q, const Letter* s, int score_cutoff, bool left, uint32_t match_mask, int sid, bool chunked, unsigned hamming_filter_id) {
	int shift = 0;
	while (mask != 0) {
		int i = ctz(mask);
		if (verify_hit(q + i + shift, s + i + shift, score_cutoff, left, match_mask >> (i + shift), sid, chunked, hamming_filter_id))
			return true;
		mask >>= i + 1;
		shift += i + 1;
	}
	return false;
}

static inline bool left_most_filter(const Sequence &query,
	const Letter* subject,
	const int seed_offset,
	const int seed_len,
	const Context &context,
	bool first_shape,
	const int shape_id,
	int score_cutoff,
	bool chunked,
	unsigned hamming_filter_id)
{
	constexpr int WINDOW_LEFT = 16, WINDOW_RIGHT = 32;

	Loc d = std::max(seed_offset - WINDOW_LEFT, 0), window_left = std::min(WINDOW_LEFT, seed_offset);
	const Letter *q = query.data() + d, *s = subject + d;
	Loc window = query.length() - d;
	window = std::min(window, window_left + 1 + WINDOW_RIGHT);

	const Sequence subject_clipped = Util::Seq::clip(s, window, window_left);
	window -= Loc(s + window - subject_clipped.end());

	d = Loc(subject_clipped.data() - s);
	q += d;
	s += d;
	window_left -= d;
	window -= d;

	const uint64_t match_mask = reduced_match(q, s, window),
		query_seed_mask = ~seed_mask(q, window);

	const uint32_t len_left = window_left + seed_len - 1,
		match_mask_left = ((1llu << len_left) - 1) & match_mask,
		query_mask_left = ((1llu << len_left) - 1) & query_seed_mask;

	const uint32_t left_hit = context.current_matcher.hit(match_mask_left, len_left) & query_mask_left;

	if (first_shape && !chunked)
		return left_hit == 0 || !verify_hits(left_hit, q, s, score_cutoff, true, match_mask_left, shape_id, chunked, hamming_filter_id);

	const uint32_t len_right = window - window_left - 1,
		match_mask_right = uint32_t(match_mask >> (window_left + 1)),
		query_mask_right = uint32_t(query_seed_mask >> (window_left + 1));
	
	const PatternMatcher& right_matcher = chunked ? context.current_matcher : context.previous_matcher;
	const uint32_t right_hit = right_matcher.hit(match_mask_right, len_right) & query_mask_right;

	return (left_hit == 0 || !verify_hits(left_hit, q, s, score_cutoff, true, match_mask_left, shape_id, chunked, hamming_filter_id))
		&& (right_hit == 0 || !verify_hits(right_hit, q + window_left + 1, s + window_left + 1, score_cutoff, false, match_mask_right, shape_id, chunked, hamming_filter_id));
}

}
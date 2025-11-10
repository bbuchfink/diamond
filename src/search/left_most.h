/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include "search.h"
#include "sse_dist.h"
#include "util/sequence/sequence.h"
#include "basic/shape_config.h"
#include "data/seed_histogram.h"
#include "hamming/finger_print.h"

namespace Search { namespace DISPATCH_ARCH {

static inline bool verify_hit(const Letter* q, const Letter* s, int score_cutoff, bool left, uint32_t match_mask, int sid, bool chunked, unsigned hamming_filter_id, PackedSeed seedp_mask) {
	using ::DISPATCH_ARCH::FingerPrint;
	if (chunked) {
		if ((shapes[sid].mask_ & match_mask) == shapes[sid].mask_) {
			PackedSeed seed;
			if (!shapes[sid].set_seed(seed, s))
				return false;
			if (left && !current_range.lower_or_equal(seed_partition(seed, seedp_mask)))
				return false;
			if (!left && !current_range.lower(seed_partition(seed, seedp_mask)))
				return false;
		}
	}
	array<char, 48> fq, fs;
	FingerPrint::load(q, &fq);
	FingerPrint::load(s, &fs);
	const unsigned id = FingerPrint(fq).match(FingerPrint(fs));
	return id >= hamming_filter_id;
}

static inline bool verify_hits(uint32_t mask, const Letter* q, const Letter* s, int score_cutoff, bool left, uint32_t match_mask, int sid, bool chunked, unsigned hamming_filter_id, PackedSeed seedp_mask) {
	int shift = 0;
	while (mask != 0) {
		int i = ctz(mask);
		if (verify_hit(q + i + shift, s + i + shift, score_cutoff, left, match_mask >> (i + shift), sid, chunked, hamming_filter_id, seedp_mask))
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
		return left_hit == 0 || !verify_hits(left_hit, q, s, score_cutoff, true, match_mask_left, shape_id, chunked, hamming_filter_id, context.seedp_mask);

	const uint32_t len_right = window - window_left - 1,
		match_mask_right = uint32_t(match_mask >> (window_left + 1)),
		query_mask_right = uint32_t(query_seed_mask >> (window_left + 1));
	
	const PatternMatcher& right_matcher = chunked ? context.current_matcher : context.previous_matcher;
	const uint32_t right_hit = right_matcher.hit(match_mask_right, len_right) & query_mask_right;

	return (left_hit == 0 || !verify_hits(left_hit, q, s, score_cutoff, true, match_mask_left, shape_id, chunked, hamming_filter_id, context.seedp_mask))
		&& (right_hit == 0 || !verify_hits(right_hit, q + window_left + 1, s + window_left + 1, score_cutoff, false, match_mask_right, shape_id, chunked, hamming_filter_id, context.seedp_mask));
}

}}
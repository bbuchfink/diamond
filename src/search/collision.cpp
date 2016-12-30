/****
Copyright (c) 2016, Benjamin Buchfink
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/

#include "collision.h"
#include "align_range.h"

// #define NO_COLLISION_FILTER

bool verify_hit(const Letter *query, const Letter *subject, unsigned sid)
{
	const Finger_print fq(query), fs(subject);
	if (fq.match(fs) < config.min_identities)
		return false;
	unsigned delta, len;
	return stage2_ungapped(query, subject, sid, delta, len) >= config.min_ungapped_raw_score;
}

inline bool match_shape_mask(const uint64_t mask, const uint64_t shape_mask)
{
	return (mask & shape_mask) == shape_mask;
}

inline bool is_lower_chunk(const Letter *subject, unsigned sid)
{
	Packed_seed seed;
	shapes[sid].set_seed(seed, subject);
	return current_range.lower(seed_partition(seed));
}

inline bool is_lower_or_equal_chunk(const Letter *subject, unsigned sid)
{
	Packed_seed seed;
	shapes[sid].set_seed(seed, subject);
	return current_range.lower_or_equal(seed_partition(seed));
}

inline bool is_high_frequency(const Letter *subject, unsigned sid, bool previous_shape)
{
	return frequent_seeds.get(subject, sid);
}

inline bool shape_collision_right(uint64_t mask, uint64_t shape_mask, const Letter *subject, unsigned sid)
{
	if (!match_shape_mask(mask, shape_mask)) return false;
	return is_lower_chunk(subject, sid)
		&& !is_high_frequency(subject, sid, false);
}

inline bool shape_collision_left(uint64_t mask, uint64_t shape_mask, const Letter *subject, unsigned sid, bool chunked)
{
	if (!match_shape_mask(mask, shape_mask)) return false;
	return (!chunked || is_lower_or_equal_chunk(subject, sid))
		&& !is_high_frequency(subject, sid, false);
}

inline bool previous_shape_collision(uint64_t mask, uint64_t shape_mask, const Letter *subject, unsigned sid)
{
	if (!match_shape_mask(mask, shape_mask)) return false;
	return !is_high_frequency(subject, sid, true);
}

bool is_primary_hit(const Letter *query,
	const Letter *subject,
	const unsigned seed_offset,
	const unsigned sid,
	const unsigned len)
{
#ifdef NO_COLLISION_FILTER
	return true;
#endif
	assert(len > 0 && len <= config.window * 2);
	const bool chunked(config.lowmem > 1);
	uint64_t mask = reduced_match32(query, subject, len);
	unsigned i = 0;
	uint64_t current_mask = shapes[sid].mask_;
	unsigned shape_len = len - shapes[0].length_ + 1;
	while (i < shape_len) {
		if (len - i > 32)
			mask |= reduced_match32(query + 32, subject + 32, len - i - 32) << 32;
		for (unsigned j = 0; j < 32 && i < shape_len; ++j) {
			/*cout << sequence(&query[j], len - j) << endl;
			cout << sequence(&subject[j], len - j) << endl;
			print_binary(mask);
			cout << endl;*/
			for (unsigned k = 0; k < sid; ++k)
				if (previous_shape_collision(mask, shapes[k].mask_, &subject[j], k) && verify_hit(&query[j], &subject[j], k)) {
					//cout << "k=" << k << endl;
					return false;
				}
			if (i < seed_offset && shape_collision_left(mask, current_mask, &subject[j], sid, chunked) && verify_hit(&query[j], &subject[j], sid))
				return false;
			if (chunked && i > seed_offset && shape_collision_right(mask, current_mask, &subject[j], sid) && verify_hit(&query[j], &subject[j], sid))
				return false;
			++i;
			mask >>= 1;
		}
		query += 32;
		subject += 32;
	}
	return true;
}
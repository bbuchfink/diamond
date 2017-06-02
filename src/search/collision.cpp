/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
	if (config.algo == Config::double_indexed)
		shapes[sid].set_seed(seed, subject);
	else
		shapes[sid].set_seed_shifted(seed, subject);
	return current_range.lower(seed_partition(seed));
}

inline bool is_lower_or_equal_chunk(const Letter *subject, unsigned sid)
{
	Packed_seed seed;
	if (config.algo == Config::double_indexed)
		shapes[sid].set_seed(seed, subject);
	else
		shapes[sid].set_seed_shifted(seed, subject);
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
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
#include <assert.h>
#include <stdint.h>
#include <limits.h>
#include <bitset>
#include <algorithm>
#include "../intrin.h"

struct PatternMatcher {

	PatternMatcher(const uint32_t* begin, const uint32_t* end) :
		min_len_(32)
	{
		const size_t count = end - begin;
		uint32_t max_len = 0;
		for (size_t i = 0; i < count; ++i) {
			assert(begin[i] != 0);
			uint32_t len = 32 - clz(begin[i]);
			max_len = std::max(max_len, len);
			min_len_ = std::min(min_len_, len);
		}
		suffix_mask_ = (1 << max_len) - 1;
		std::fill(table_, table_ + SIZE, '\0');
		for (uint32_t s = 0; s <= suffix_mask_; ++s) {
			for (size_t i = 0; i < count; ++i)
				if ((s & begin[i]) == begin[i])
					//table_[s] = true;
					table_[s] = 1;
		}
	}

	uint32_t hit(uint32_t h, uint32_t len) const {
		if (len < min_len_)
			return 0;
		const uint32_t mask = suffix_mask_, end = len - min_len_ + 1;
		uint32_t r = 0;
		for (uint32_t i = 0; i < end; ++i) {
			r |= uint32_t(table_[h & mask]) << i;
			h >>= 1;
		}
		return r;
	}

private:
	
	static constexpr uint64_t SIZE = 1 << MAX_SHAPE_LEN;

	uint32_t min_len_, suffix_mask_;
	//std::bitset<SIZE> table_;
	unsigned char table_[SIZE];
	
};

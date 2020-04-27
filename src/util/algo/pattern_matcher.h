#ifndef PATTERN_MATCHER_H_
#define PATTERN_MATCHER_H_

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

#endif
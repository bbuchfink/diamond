#ifndef PATTERN_MATCHER_H_
#define PATTERN_MATCHER_H_

#include <stdint.h>
#include <bitset>
#include <vector>
#include <algorithm>
#include "../intrin.h"

struct PatternMatcher {

	PatternMatcher(const uint64_t* begin, const uint64_t* end) :
		count_(end - begin),
		max_len_(0)
	{
		for (size_t i = 0; i < count_; ++i) {
			len_.push_back(64 - clz(begin[i]));
			max_len_ = std::max(max_len_, len_.back());
		}
		suffix_count_ = 1 << max_len_;
		for (uint64_t s = 0; s < suffix_count_; ++s) {
			for (size_t i = 0; i < count_; ++i)
				if (hit(s, begin[i], len_[i]))
					table_[s] = true;
		}
	}

	uint32_t hit(uint32_t h) const {

	}

private:

	bool hit(uint64_t suffix, uint64_t pattern, size_t len) {
		for (size_t i = 0; i < max_len_ - len; ++i) {
			if ((suffix & pattern) == pattern)
				return true;
			suffix >>= 1;
		}
		return false;
	}

	constexpr uint64_t SIZE = 1 << MAX_SHAPE_LEN;

	const size_t count_;
	size_t max_len_;
	uint64_t suffix_count_;
	std::vector<size_t> len_;
	std::bitset<SIZE> table_;
	
};

#endif
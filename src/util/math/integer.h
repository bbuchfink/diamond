#include <limits.h>
#include <stdint.h>
#include <algorithm>
#include <math.h>

#ifndef INTEGER_H_
#define INTEGER_H_

static inline int8_t saturated_add(int8_t x, int8_t y) {
	return (int8_t)std::max(int32_t(x) + int32_t(y), (int32_t)std::numeric_limits<int8_t>::min());
}

static inline int16_t saturated_add(int16_t x, int16_t y) {
	return (int16_t)std::max(int32_t(x) + int32_t(y), (int32_t)std::numeric_limits<int16_t>::min());
}

static inline int32_t saturated_add(int32_t x, int32_t y) {
	return x + y;
}

static inline size_t bit_length(size_t x) {
	return ceil(log(x) / log(2)) + 1;
}

static inline uint64_t next_power_of_2(double x)
{
	return 1llu << uint64_t(ceil(log(x) / log(2)));
}

#endif
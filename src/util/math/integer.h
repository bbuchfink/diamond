#include <limits.h>
#include <stdint.h>
#include <algorithm>

#ifndef INTEGER_H_
#define INTEGER_H_

inline int16_t saturated_add(int16_t x, int16_t y) {
	return (int16_t)std::max(int32_t(x) + int32_t(y), (int32_t)std::numeric_limits<int16_t>::min());
}

inline int32_t saturated_add(int32_t x, int32_t y) {
	return x + y;
}

#endif
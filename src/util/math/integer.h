/****
Copyright © 2012-2026 Benjamin J. Buchfink <buchfink@gmail.com>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

	http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
****/
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <cmath>
#include <math.h>
#include <limits.h>
#include <stdint.h>
#include <assert.h>
#include <limits>
#include "../intrin.h"

static inline int8_t saturated_add(int8_t x, int8_t y) {
	return (int8_t)std::max(int32_t(x) + int32_t(y), (int32_t)std::numeric_limits<int8_t>::min());
}

static inline int16_t saturated_add(int16_t x, int16_t y) {
	return (int16_t)std::max(int32_t(x) + int32_t(y), (int32_t)std::numeric_limits<int16_t>::min());
}

static inline int32_t saturated_add(int32_t x, int32_t y) {
	return x + y;
}

static inline int bit_length(uint64_t x) {
	assert(x > 0);
	return 64 - clz(x);
}

inline size_t next_pow2(size_t x) noexcept {
	if (x <= 1) return 1;
	if (x > (std::numeric_limits<size_t>::max() >> 1)) return 0;
#if SIZE_MAX == 0xFFFFFFFFu
	unsigned lz = clz(static_cast<uint32_t>(x - 1));
	const unsigned W = 32;
#elif SIZE_MAX == 0xFFFFFFFFFFFFFFFFull
	unsigned lz = clz(static_cast<uint64_t>(x - 1));
	const unsigned W = 64;
#else
#error "Unsupported size_t width"
#endif
	return size_t{ 1 } << (W - lz);
}

inline size_t next_pow2(double x) noexcept {
	return next_pow2(static_cast<size_t>(std::ceil(x)));
}

template<typename It, int N>
bool next_combination(It begin, It end) {
	It i = begin;
	while (i < end) {
		if (*i < N - 1) {
			++(*i);
			return true;
		}
		else {
			*i = 0;
			++i;
		}
	}
	return false;
}

template<typename I>
I power(I x, I p)
{
	if (p == 0) return 1;
	if (p == 1) return x;

	const I t = power(x, p / 2);
	if (p % 2 == 0)
		return t * t;
	else
		return x * t * t;
}

template<typename Int>
int digits(Int x, int base) {
	int d = 0;
	while (x > 0) {
		x /= base;
		++d;
	}
	return d;
}
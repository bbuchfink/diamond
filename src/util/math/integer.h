/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once
#include <cmath>
#include <math.h>
#include <limits.h>
#include <stdexcept>
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
	return x == 0 ? 0 : 64 - clz(x);
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

inline uint64_t gb_to_bytes(double gb) {
	long double bytes = static_cast<long double>(gb) * 1000000000.0L;
	if (!std::isfinite(gb) ||
		bytes < 0 ||
		bytes > static_cast<long double>(std::numeric_limits<uint64_t>::max())) {
		throw std::overflow_error("block size does not fit in uint64_t");
	}
	return static_cast<uint64_t>(bytes);
}
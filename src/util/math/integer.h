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

size_t next_pow2(size_t x) noexcept {
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

size_t next_pow2(double x) noexcept {
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
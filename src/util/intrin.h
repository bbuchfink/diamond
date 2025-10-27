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
#include <stdint.h>
#ifdef _MSC_VER
#include <intrin.h>
#endif

static inline unsigned popcount32(unsigned x)
{
#ifdef _MSC_VER
	return __popcnt(x);
#else
	return __builtin_popcount(x);
#endif
}

static inline unsigned popcount64(unsigned long long x)
{
#ifdef _MSC_VER
	return (unsigned)__popcnt64(x);
#else
	return __builtin_popcountll(x);
#endif
}

static inline int ctz(uint32_t x)
{
#ifdef _MSC_VER
	unsigned long i;
	unsigned char c = _BitScanForward(&i, x);
	return i;
#else
	return __builtin_ctz(x);
#endif
}


static inline int ctz(uint64_t x)
{
#ifdef _MSC_VER
	if (x)
		return (int)__popcnt64((x ^ (x - 1)) >> 1);
	else
		return CHAR_BIT * sizeof(x);
#else
	return __builtin_ctzll(x);
#endif
}

static inline int clz(uint64_t x) {
#ifdef _MSC_VER
	unsigned long i;
	unsigned char c = _BitScanReverse64(&i, x);
	return 63 - (int)i;
#else
	return __builtin_clzll(x);
#endif
}

static inline int clz(uint32_t x) {
#ifdef _MSC_VER
	unsigned long i;
	unsigned char c = _BitScanReverse(&i, x);
	return 31 - (int)i;
#else
	return __builtin_clz(x);
#endif
}
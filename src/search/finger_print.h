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
#include "util/simd.h"
#include "basic/value.h"
#include <cstring>

#ifdef __SSE2__

struct Byte_finger_print_48
{
	Byte_finger_print_48(const Letter *q) :
#ifdef SEQ_MASK
		r1(letter_mask(_mm_loadu_si128((__m128i const*)(q - 16)))),
		r2(letter_mask(_mm_loadu_si128((__m128i const*)(q)))),
		r3(letter_mask(_mm_loadu_si128((__m128i const*)(q + 16))))
#else
		r1(_mm_loadu_si128((__m128i const*)(q - 16))),
		r2(_mm_loadu_si128((__m128i const*)(q))),
		r3(_mm_loadu_si128((__m128i const*)(q + 16)))
#endif
	{}
	static uint64_t match_block(__m128i x, __m128i y)
	{
		return (uint64_t)_mm_movemask_epi8(_mm_cmpeq_epi8(x, y));
	}
	unsigned match(const Byte_finger_print_48 &rhs) const
	{
		return popcount64(match_block(r3, rhs.r3) << 32 | match_block(r1, rhs.r1) << 16 | match_block(r2, rhs.r2));
	}
	alignas(16) __m128i r1, r2, r3;
};

#elif defined(__ARM_NEON)

struct Byte_finger_print_48
{
	Byte_finger_print_48(const Letter *q) :
#ifdef SEQ_MASK
		r1(letter_mask(vld1q_s8(q - 16))),
		r2(letter_mask(vld1q_s8(q))),
		r3(letter_mask(vld1q_s8(q + 16)))
#else
		r1(vld1q_s8(q - 16)),
		r2(vld1q_s8(q)),
		r3(vld1q_s8(q + 16)),
#endif
	{}
	unsigned match(const Byte_finger_print_48 &rhs) const
	{
		const uint8x16_t ONES = vdupq_n_u8(1);
		uint8x16_t s1 = vandq_u8(vceqq_s8(r1, rhs.r1), ONES);
		uint8x16_t s2 = vandq_u8(vceqq_s8(r2, rhs.r2), ONES);
		uint8x16_t s3 = vandq_u8(vceqq_s8(r3, rhs.r3), ONES);
		uint16x8_t acc = vdupq_n_u16(0);
		acc = vpadalq_u8(acc, s1);
		acc = vpadalq_u8(acc, s2);
		acc = vpadalq_u8(acc, s3);
		return vhsumq_u16(acc);
	}
	alignas(16) int8x16_t r1, r2, r3;
};

#else

struct Byte_finger_print_48
{
	Byte_finger_print_48()
	{}
	Byte_finger_print_48(const Letter *q)
	{
		memcpy(r, q - 16, 48);
#ifdef SEQ_MASK
		for (int i = 0; i < 48; ++i)
			r[i] &= LETTER_MASK;
#endif
	}
	unsigned match(const Byte_finger_print_48 &rhs) const
	{
		unsigned n = 0;
		for (unsigned i = 0; i < 48; ++i)
			if (r[i] == rhs.r[i])
				++n;
		return n;
	}
	Letter r[48];
};

#endif

typedef Byte_finger_print_48 FingerPrint;
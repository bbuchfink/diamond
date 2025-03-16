/****
DIAMOND protein aligner
Copyright (C) 2016-2024 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
						
Code developed by Benjamin Buchfink <buchfink@gmail.com>
Arm NEON port contributed by Martin Larralde <martin.larralde@embl.de>

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
#include "util/simd.h"
#include "basic/value.h"
#include <cstring>

#ifdef __AVX2__

struct Byte_finger_print_32
{
	Byte_finger_print_32(const Letter* q) :
#ifdef SEQ_MASK
		r1(letter_mask(_mm256_loadu_si256((__m256i const*)(q - 16))))
#else
		r1(_mm256_loadu_si256((__m256i const*)(q - 16)))
#endif
	{}
	static uint64_t match_block(__m256i x, __m256i y)
	{
		return (uint64_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(x, y));
	}
	unsigned match(const Byte_finger_print_32& rhs) const
	{
		return popcount64(match_block(r1, rhs.r1));
	}
	__m256i r1;
};

struct Byte_finger_print_64
{
	Byte_finger_print_64(const Letter* q) :
#ifdef SEQ_MASK
		r1(letter_mask(_mm256_loadu_si256((__m256i const*)(q - 32)))),
		r2(letter_mask(_mm256_loadu_si256((__m256i const*)q)))
#else
		r1(_mm256_loadu_si256((__m256i const*)(q - 32))),
		r2(_mm256_loadu_si256((__m256i const*)q))
#endif
	{}
	static uint64_t match_block(__m256i x, __m256i y)
	{
		return (uint64_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(x, y));
	}
	unsigned match(const Byte_finger_print_64& rhs) const
	{
		return popcount64(match_block(r1, rhs.r1) << 32 | match_block(r2, rhs.r2));
	}
	__m256i r1, r2;
};

#endif

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
		//printf("%llx\n", q);
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
	//char r[32];
};

#endif

typedef Byte_finger_print_48 FingerPrint;
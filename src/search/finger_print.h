/****
DIAMOND protein aligner
Copyright (C) 2016-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
						
Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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
#include "../util/simd.h"
#include "../basic/config.h"

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
	bool operator==(const Byte_finger_print_48& rhs) const {
		return match(rhs) >= config.min_identities;
	}
	alignas(16) __m128i r1, r2, r3;
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
	bool operator==(const Byte_finger_print_48& rhs) const {
		return match(rhs) >= config.min_identities;
	}
	Letter r[48];
	//char r[32];
};

#endif

#ifdef __AVX2__
typedef Byte_finger_print_48 Finger_print;
#else
typedef Byte_finger_print_48 Finger_print;
#endif
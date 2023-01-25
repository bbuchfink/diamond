/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen
						
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
#include <assert.h>
#include "../basic/reduction.h"
#include "../basic/value.h"
#include "../util/simd.h"

#ifdef __SSSE3__
static inline __m128i reduce_seq_ssse3(const Letter *seq, const Letter* map)
{
	const __m128i *row = reinterpret_cast<const __m128i*>(map);
	__m128i s = _mm_loadu_si128((const __m128i*)seq);
#ifdef SEQ_MASK
	s = letter_mask(s);
#endif
	__m128i high_mask = _mm_slli_epi16(_mm_and_si128(s, _mm_set1_epi8('\x10')), 3);
	__m128i seq_low = _mm_or_si128(s, high_mask);
	__m128i seq_high = _mm_or_si128(s, _mm_xor_si128(high_mask, _mm_set1_epi8('\x80')));

	__m128i r1 = _mm_loadu_si128(row);
	__m128i r2 = _mm_loadu_si128(row+1);
	__m128i s1 = _mm_shuffle_epi8(r1, seq_low);
	__m128i s2 = _mm_shuffle_epi8(r2, seq_high);
	return _mm_or_si128(s1, s2);
}
#endif

#ifdef __SSE2__
static inline __m128i reduce_seq_generic(const Letter *seq, const Letter *map)
{
	alignas(16) uint8_t d[16];
	for (unsigned i = 0; i < 16; ++i)
		d[i] = map[(long)letter_mask(*(seq++))];
	return _mm_load_si128((const __m128i*)d);
}

static inline __m128i reduce_seq(const Letter *seq, const Letter *map)
{
#ifdef __SSSE3__
	return reduce_seq_ssse3(seq, map);
#else
	return reduce_seq_generic(seq, map);
#endif
}
#endif

static inline unsigned match_block_reduced(const Letter *x, const Letter *y)
{
#ifdef __SSE2__
	const __m128i r1 = reduce_seq(x, Reduction::reduction.map8());
	const __m128i r2 = reduce_seq(y, Reduction::reduction.map8b());
	return _mm_movemask_epi8(_mm_cmpeq_epi8(r1, r2));
#else
	unsigned r = 0;
	for (int i = 15; i >= 0; --i) {
		r <<= 1;
		const Letter lx = letter_mask(x[i]), ly = letter_mask(y[i]);
		if (!is_amino_acid(lx) || !is_amino_acid(ly))
			continue;
		if (Reduction::reduction(lx) == Reduction::reduction(ly))
			r |= 1;
	}
	return r;
#endif
}

static inline uint64_t reduced_match32(const Letter* q, const Letter* s, unsigned len)
{
	uint64_t x = match_block_reduced(q + 16, s + 16) << 16 | match_block_reduced(q, s);
	if (len < 32)
		x &= (1 << len) - 1;
	return x;
}

static inline uint64_t reduced_match(const Letter* q, const Letter* s, int len) {
	assert(len <= 64);
	if (len < 64) {
		const uint64_t mask = (1llu << len) - 1;
		uint64_t m = match_block_reduced(q, s);
		if (len <= 16)
			return m & mask;
		m |= (uint64_t)match_block_reduced(q + 16, s + 16) << 16;
		if (len <= 32)
			return m & mask;
		m |= (uint64_t)match_block_reduced(q + 32, s + 32) << 32;
		if (len <= 48)
			return m & mask;
		m |= (uint64_t)match_block_reduced(q + 48, s + 48) << 48;
		return m & mask;
	}
	else
		return (uint64_t)match_block_reduced(q, s) | ((uint64_t)match_block_reduced(q + 16, s + 16) << 16) | ((uint64_t)match_block_reduced(q + 32, s + 32) << 32) | ((uint64_t)match_block_reduced(q + 48, s + 48) << 48);
}

static inline uint64_t seed_mask(const Letter* s, int len) {
#ifdef __SSE2__
	assert(len <= 64);
	uint64_t mask = 0;
	const __m128i m = _mm_set1_epi8(SEED_MASK);
	for (int i = 0; i < len; i += 16) {
		const __m128i mask_bits = _mm_and_si128(_mm_loadu_si128((const __m128i*)s), m);
		mask |= (uint64_t)_mm_movemask_epi8(mask_bits) << i;
		s += 16;
	}
	if (len < 64)
		mask &= (1llu << len) - 1;
	return mask;
#else
	uint64_t mask = 0;
	for (int i = 0; i < len; ++i) {
		if ((s[i] & SEED_MASK) != 0)
			mask |= 1llu << i;
	}
	return mask;
#endif
}

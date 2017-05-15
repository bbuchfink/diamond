/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#ifndef SSE_DIST_H_
#define SSE_DIST_H_

#include "../basic/reduction.h"
#include "../basic/value.h"
#include "../util/simd.h"

#ifdef __SSSE3__
inline __m128i reduce_seq_ssse3(const __m128i &seq)
{
	const __m128i *row = reinterpret_cast<const __m128i*>(Reduction::reduction.map8());
	__m128i high_mask = _mm_slli_epi16(_mm_and_si128(seq, _mm_set1_epi8('\x10')), 3);
	__m128i seq_low = _mm_or_si128(seq, high_mask);
	__m128i seq_high = _mm_or_si128(seq, _mm_xor_si128(high_mask, _mm_set1_epi8('\x80')));

	__m128i r1 = _mm_load_si128(row);
	__m128i r2 = _mm_load_si128(row+1);
	__m128i s1 = _mm_shuffle_epi8(r1, seq_low);
	__m128i s2 = _mm_shuffle_epi8(r2, seq_high);
	return _mm_or_si128(s1, s2);
}
#endif

#ifdef __SSE2__
inline __m128i reduce_seq_generic(const __m128i &seq)
{
	__m128i r;
	Letter* s = (Letter*)&seq;
	uint8_t* d = (uint8_t*)&r;
	for(unsigned i=0;i<16;++i)
		*(d++) = Reduction::reduction(*(s++));
	return r;
}

inline __m128i reduce_seq(const __m128i &seq)
{
#ifdef __SSSE3__
	return reduce_seq_ssse3(seq);
#else
	return reduce_seq_generic(seq);
#endif
}
#endif

inline unsigned match_block_reduced(const Letter *x, const Letter *y)
{
#ifdef __SSE2__
	static const __m128i mask = _mm_set1_epi8(0x7F);
	__m128i r1 = _mm_loadu_si128 ((__m128i const*)(x));
	__m128i r2 = _mm_loadu_si128 ((__m128i const*)(y));
	r2 = _mm_and_si128(r2, mask);
	r1 = reduce_seq(r1);
	r2 = reduce_seq(r2);
	return _mm_movemask_epi8(_mm_cmpeq_epi8(r1, r2));
#else
	unsigned r = 0;
	for (int i = 15; i >= 0; --i) {
		r <<= 1;
		if (Reduction::reduction(x[i]) == Reduction::reduction(y[i]))
			r |= 1;		
	}
	return r;
#endif
}

inline uint64_t reduced_match32(const Letter* q, const Letter *s, unsigned len)
{
	uint64_t x = match_block_reduced(q + 16, s + 16) << 16 | match_block_reduced(q, s);
	if(len < 32)
		x &= (1 << len) - 1;
	return x;
}

#endif /* SSE_DIST_H_ */

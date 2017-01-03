/****
Copyright (c) 2013-2016, Benjamin Buchfink
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
	if(config.have_ssse3) {
#ifdef __SSSE3__
		return reduce_seq_ssse3(seq);
#else
		return reduce_seq_generic(seq);
#endif
	} else
		return reduce_seq_generic(seq);
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

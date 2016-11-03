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

inline unsigned match_block(const Letter *x, const Letter *y)
{
	static const __m128i mask = _mm_set1_epi8(0x7F);
	__m128i r1 = _mm_loadu_si128 ((__m128i const*)(x));
	__m128i r2 = _mm_loadu_si128 ((__m128i const*)(y));
	r2 = _mm_and_si128(r2, mask);
	return _mm_movemask_epi8(_mm_cmpeq_epi8(r1, r2));
}

inline unsigned popcount32(unsigned x)
{
#ifdef _MSC_VER
	return __popcnt(x);
	//return (unsigned)popcount_3(x);
#else
	return __builtin_popcount(x);
#endif
}

inline unsigned popcount64(unsigned long long x)
{
#ifdef _MSC_VER
	return (unsigned)__popcnt64(x);
#else
	return __builtin_popcountll(x);
#endif
}

inline unsigned fast_match(const Letter *q, const Letter *s)
{
	return popcount32(match_block(q - 8, s - 8) << 16 | match_block(q + 8, s + 8));
}

struct Masked {};

struct Int_finger_print
{
	Int_finger_print(const Letter *q) :
		r1(*(uint64_t*)(q-8)),
		r2(*(uint64_t*)(q)),
		r3(*(uint64_t*)(q + 8)),
		r4(*(uint64_t*)(q +16))
	{ }
	Int_finger_print(const Letter *q, Masked) :
		r1(*(uint64_t*)(q - 8)),
		r2(*(uint64_t*)(q)),
		r3(*(uint64_t*)(q + 8)),
		r4(*(uint64_t*)(q + 16))
	{ }
	unsigned match(const Int_finger_print &rhs) const
	{
		return popcount64(r1 & rhs.r1) + popcount64(r2 & rhs.r2) + popcount64(r3 & rhs.r3) + popcount64(r4 & rhs.r4);
	}
	uint64_t r1, r2,r3,r4;
};

#ifdef __SSE2__

struct Byte_finger_print
{
	Byte_finger_print(const Letter *q) :
		r1(_mm_loadu_si128((__m128i const*)(q - 8))),
		r2(_mm_loadu_si128((__m128i const*)(q + 8)))
	{ }
	Byte_finger_print(const Letter *q, Masked):
		r1 (_mm_and_si128(_mm_loadu_si128((__m128i const*)(q - 8)), _mm_set1_epi8('\x7f'))),
		r2 (_mm_and_si128(_mm_loadu_si128((__m128i const*)(q + 8)), _mm_set1_epi8('\x7f')))
	{ }
	static unsigned match_block(__m128i x, __m128i y)
	{
		return _mm_movemask_epi8(_mm_cmpeq_epi8(x, y));
	}
	unsigned match(const Byte_finger_print &rhs) const
	{
		//return popcount_3(match_block(r1, rhs.r1) << 16 | match_block(r2, rhs.r2));
		return popcount32(match_block(r1, rhs.r1) << 16 | match_block(r2, rhs.r2));
	}
	__m128i r1, r2;
};

#else

struct Byte_finger_print
{
	Byte_finger_print()
	{}
	Byte_finger_print(const Letter *q)
	{
		//printf("%llx\n", q);
		memcpy(r, q - 8, 32);
	}
	unsigned match(const Byte_finger_print &rhs) const
	{
		unsigned n = 0;
		for (unsigned i = 0; i < 32; ++i)
			if (r[i] == rhs.r[i])
				++n;
		return n;
	}
	Letter r[32];
	//char r[32];
};

#endif

struct Halfbyte_finger_print_naive
{
	Halfbyte_finger_print_naive(const Letter *q) :
		r1(reduce(q-8)),
		r2(reduce(q+8))
	{ }
	static unsigned match_block(uint64_t x, uint64_t y)
	{
		uint64_t v = ~(x ^ y);
		v &= v >> 1;
		v &= 0x5555555555555555LL;
		v &= v >> 2;
		v &= 0x1111111111111111LL;
		return (unsigned)popcount64(v);
	}
	unsigned match(const Halfbyte_finger_print_naive &rhs) const
	{
		return match_block(r1, rhs.r1) + match_block(r2, rhs.r2);
	}
	static uint64_t reduce(const Letter *q)
	{
		uint64_t x;
		for (unsigned i = 0; i < 16; ++i)
			x = (x << 4) | reduction(q[i]);
		return x;
	}
	static const Reduction reduction;
	uint64_t r1, r2;
};

struct Halfbyte_finger_print
{
	Halfbyte_finger_print(const Letter *q) :
		r(_mm_set_epi32(reduce(q - 8), reduce(q), reduce(q + 8), reduce(q+16)))
	{ }
	static unsigned match_block(uint64_t x, uint64_t y)
	{
		uint64_t v = ~(x ^ y);
		v &= v >> 1;
		v &= 0x5555555555555555LL;
		v &= v >> 2;
		v &= 0x1111111111111111LL;
		return (unsigned)popcount64(v);
	}
	static unsigned get_mask(__m128i x, __m128i y)
	{
		return _mm_movemask_epi8(_mm_cmpeq_epi8(x, y));
	}
	unsigned match(const Halfbyte_finger_print &rhs) const
	{
		return popcount32(get_mask(_mm_and_si128(r,_mm_set1_epi8('\xf0')), _mm_and_si128(rhs.r, _mm_set1_epi8('\xf0')))<<16
			| get_mask(_mm_and_si128(r, _mm_set1_epi8('\x0f')), _mm_and_si128(rhs.r, _mm_set1_epi8('\x0f'))));
	}
	static int reduce(const Letter *q)
	{
		int x = 0;
		for (unsigned i = 0; i < 8; ++i)
			x = (x << 4) | reduction(q[i]);
		return x;
	}
	static const Reduction reduction;
	__m128i r;
};

typedef Byte_finger_print Finger_print;

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
	uint64_t x = match_block_reduced(q+16, s+16)<<16 | match_block_reduced(q,s);
	if(len < 32)
		x &= (1 << len) - 1;
	return x;
}

#endif /* SSE_DIST_H_ */

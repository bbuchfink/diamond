/****
DIAMOND protein aligner
Copyright (C) 2013-2019 Benjamin Buchfink <buchfink@gmail.com>

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

#ifndef SCORE_VECTOR_H_
#define SCORE_VECTOR_H_

#include <algorithm>
#include <limits.h>
#include "../util/simd.h"
#include "../basic/score_matrix.h"
#include "../util/simd/vector.h"

template<typename _score>
struct score_traits
{
	static const unsigned channels = 1;
	enum { zero = 0, byte_size = 4 };
	typedef bool Mask;
};

template<>
struct score_traits<uint8_t>
{
	//static const unsigned channels = 16;
	enum { channels = 16, zero = 0x00, byte_size = 1 };
	typedef uint16_t Mask;
};

namespace DISPATCH_ARCH {

template<typename _t>
struct ScoreTraits
{
	typedef void Score;
	enum { CHANNELS = 0 };
};


template<typename _score>
struct score_vector
{ };

}

namespace DISPATCH_ARCH {

#ifdef __SSE2__

template<>
struct score_vector<uint8_t>
{

	typedef uint8_t Score;
	enum { CHANNELS = 16 };

	score_vector()
	{
		//data_ = _mm_set1_epi8(score_traits<uint8_t>::zero);
		data_ = _mm_setzero_si128();
	}

	explicit score_vector(uint8_t x):
		data_ (::SIMD::_mm_set1_epi8((char)x))
	{ }

	explicit score_vector(__m128i data):
		data_ (data)
	{ }

	score_vector(unsigned a, uint64_t seq)
	{}

	score_vector(unsigned a, const __m128i &seq)
	{
#ifdef __SSSE3__
		set_ssse3(a, seq);
#else
		set_generic(a, seq);
#endif
	}

	score_vector(unsigned a, const __m128i &seq, const score_vector &bias)
	{
#ifdef __SSSE3__
		set_ssse3(a, seq);
#else
		set_generic(a, seq);
#endif
	}

	void set_ssse3(unsigned a, const __m128i &seq)
	{
#ifdef __SSSE3__
		const __m128i *row = reinterpret_cast<const __m128i*>(&score_matrix.matrix8u()[a << 5]);

		__m128i high_mask = _mm_slli_epi16(_mm_and_si128(seq, ::SIMD::_mm_set1_epi8('\x10')), 3);
		__m128i seq_low = _mm_or_si128(seq, high_mask);
		__m128i seq_high = _mm_or_si128(seq, _mm_xor_si128(high_mask, ::SIMD::_mm_set1_epi8('\x80')));

		__m128i r1 = _mm_loadu_si128(row);
		__m128i r2 = _mm_loadu_si128(row+1);
		__m128i s1 = _mm_shuffle_epi8(r1, seq_low);
		__m128i s2 = _mm_shuffle_epi8(r2, seq_high);
		data_ = _mm_or_si128(s1, s2);
#endif
	}
	
	void set_generic(unsigned a, const __m128i &seq)
	{
		const uint8_t* row (&score_matrix.matrix8u()[a<<5]);
		uint8_t l[16], dest[16];
		_mm_storeu_si128((__m128i*)l, seq);
		for (unsigned i = 0; i < 16; i++)
			dest[i] = row[(unsigned long)l[i]];
		data_ = _mm_loadu_si128((const __m128i*)dest);
	}

	score_vector(const uint8_t* s):
		data_ (_mm_loadu_si128(reinterpret_cast<const __m128i*>(s)))
	{ }

	score_vector operator+(const score_vector &rhs) const
	{
		return score_vector (_mm_adds_epu8(data_, rhs.data_));
	}

	score_vector operator-(const score_vector &rhs) const
	{
		return score_vector (_mm_subs_epu8(data_, rhs.data_));
	}

	score_vector& operator-=(const score_vector &rhs)
	{
		data_ = _mm_subs_epu8(data_, rhs.data_);
		return *this;
	}

	score_vector& operator &=(const score_vector& rhs) {
		data_ = _mm_and_si128(data_, rhs.data_);
		return *this;
	}

	__m128i operator==(const score_vector &rhs) const
	{
		return _mm_cmpeq_epi8(data_, rhs.data_);
	}

	void unbias(const score_vector &bias)
	{ this->operator -=(bias); }

	int operator [](unsigned i) const
	{
		uint8_t s[16];
		_mm_storeu_si128((__m128i*)s, data_);
		return (int)s[i];
	}

	void set(unsigned i, uint8_t v)
	{
		uint8_t s[16];
		_mm_storeu_si128((__m128i*)s, data_);
		s[i] = v;
		data_ = _mm_loadu_si128((const __m128i*)s);
	}

	score_vector& max(const score_vector &rhs)
	{
		data_ = _mm_max_epu8(data_, rhs.data_);
		return *this;
	}

	score_vector& min(const score_vector &rhs)
	{
		data_ = _mm_min_epu8(data_, rhs.data_);
		return *this;
	}

	friend score_vector max(const score_vector& lhs, const score_vector &rhs)
	{
		return score_vector (_mm_max_epu8(lhs.data_, rhs.data_));
	}

	friend score_vector min(const score_vector& lhs, const score_vector &rhs)
	{
		return score_vector (_mm_min_epu8(lhs.data_, rhs.data_));
	}

	uint16_t cmpeq(const score_vector &rhs) const
	{
		return _mm_movemask_epi8(_mm_cmpeq_epi8(data_, rhs.data_));
	}

	__m128i cmpeq2(const score_vector &rhs) const
	{
		return _mm_cmpeq_epi8(data_, rhs.data_);
	}

	uint16_t cmpgt(const score_vector &rhs) const
	{
		return _mm_movemask_epi8(_mm_cmpgt_epi8(data_, rhs.data_));
	}

	void store(uint8_t *ptr) const
	{
		_mm_storeu_si128((__m128i*)ptr, data_);
	}

	bool operator>(score_vector<uint8_t> cmp) const
	{
		const score_vector<uint8_t> s = *this - cmp;
#ifdef __SSE4_1__
		return _mm_testz_si128(s.data_, s.data_) == 0;
#else
		return _mm_movemask_epi8(_mm_cmpeq_epi8(s.data_, _mm_setzero_si128())) == 0xFFFF;
#endif
	}

	friend std::ostream& operator<<(std::ostream &s, score_vector v)
	{
		uint8_t x[16];
		v.store(x);
		for (unsigned i = 0; i < 16; ++i)
			printf("%3i ", (int)x[i]);
		return s;
	}

	__m128i data_;

};

#endif

template<>
struct ScoreTraits<int32_t>
{
	enum { CHANNELS = 1, BITS = 32 };
	typedef int32_t Score;
	struct TraceMask {
		static uint8_t make(int vmask, int hmask) {
			return vmask << 1 | hmask;
		}
		static uint8_t vmask(int channel) {
			return 2;
		}
		static uint8_t hmask(int channel) {
			return 1;
		}
		uint8_t gap;
		uint8_t open;
	};
	static int32_t zero()
	{
		return 0;
	}
	static int32_t zero_score()
	{
		return 0;
	}
	static int int_score(Score s)
	{
		return s;
	}
	static void saturate(int &x)
	{
		x = std::max(x, 0);
	}
	static int max_score()
	{
		return INT_MAX;
	}
	static constexpr int max_int_score() {
		return INT_MAX;
	}
};

#ifdef __SSE2__
template<>
struct ScoreTraits<score_vector<uint8_t>>
{
	enum { CHANNELS = 16 };
	static score_vector<uint8_t> zero() {
		return score_vector<uint8_t>();
	}
	static constexpr uint8_t max_score() {
		return std::numeric_limits<uint8_t>::max();
	}
};
#endif

}

static inline void store_sv(int32_t sv, int32_t *dst)
{
	*dst = sv;
}

static inline int32_t load_sv(const int32_t *x) {
	return *x;
}

#ifdef __SSE2__

template<typename _t, typename _p>
static inline void store_sv(const DISPATCH_ARCH::score_vector<_t> &sv, _p *dst)
{
#if ARCH_ID == 2
	_mm256_storeu_si256((__m256i*)dst, sv.data_);
#else
	_mm_storeu_si128((__m128i*)dst, sv.data_);
#endif
}

#endif

template<typename _sv>
static inline typename DISPATCH_ARCH::ScoreTraits<_sv>::Score extract_channel(const _sv &v, int i) {
	return v[i];
}

template<>
inline int extract_channel<int>(const int &v, int i) {
	return v;
}

template<typename _sv>
static inline void set_channel(_sv &v, int i, typename DISPATCH_ARCH::ScoreTraits<_sv>::Score x) {
	v.set(i, x);
}

template<>
inline void set_channel<int>(int &v, int i, int x) {
	v = x;
}

#endif /* SCORE_VECTOR_H_ */

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

#ifndef SCORE_VECTOR_H_
#define SCORE_VECTOR_H_

#include <limits.h>
#include "../util/simd.h"
#include "../basic/score_matrix.h"

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

template<typename _score>
struct score_vector
{ };

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

	explicit score_vector(char x):
		data_ (_mm_set(x))
	{ }

	explicit score_vector(__m128i data):
		data_ (data)
	{ }

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

		__m128i high_mask = _mm_slli_epi16(_mm_and_si128(seq, _mm_set1_epi8('\x10')), 3);
		__m128i seq_low = _mm_or_si128(seq, high_mask);
		__m128i seq_high = _mm_or_si128(seq, _mm_xor_si128(high_mask, _mm_set1_epi8('\x80')));

		__m128i r1 = _mm_load_si128(row);
		__m128i r2 = _mm_load_si128(row+1);
		__m128i s1 = _mm_shuffle_epi8(r1, seq_low);
		__m128i s2 = _mm_shuffle_epi8(r2, seq_high);
		data_ = _mm_or_si128(s1, s2);
#endif
	}
	
	void set_generic(unsigned a, const __m128i &seq)
	{
		const uint8_t* row (&score_matrix.matrix8u()[a<<5]);
		const uint8_t* seq_ptr (reinterpret_cast<const uint8_t*>(&seq));
		uint8_t* dest (reinterpret_cast<uint8_t*>(&data_));
		for(unsigned i=0;i<16;i++)
			*(dest++) = row[*(seq_ptr++)];
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

	score_vector& operator++()
	{
		data_ = _mm_adds_epu8(data_, _mm_set(1));
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
		return *(((uint8_t*)&data_)+i);
	}

	void set(unsigned i, uint8_t v)
	{
		*(((uint8_t*)&data_)+i) = v;
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

template<>
struct score_traits<int16_t>
{
	enum { channels = 8, zero = 0x8000, byte_size = 2 };
	typedef uint8_t Mask;
};

template<>
struct score_vector<int16_t>
{

	typedef int16_t Score;
	enum { CHANNELS = 8 };

	score_vector()
	{
	}

	score_vector(int x)
	{
		data_ = _mm_set1_epi16(x);
	}

	explicit score_vector(__m128i data) :
		data_(data)
	{ }

	score_vector(unsigned a, uint64_t seq)
	{
		const uint16_t* row((uint16_t*)&score_matrix.matrix16()[a << 5]);
		uint64_t b = uint64_t(row[seq & 0xff]);
		seq >>= 8;
		b |= uint64_t(row[seq & 0xff]) << 16;
		seq >>= 8;
		b |= uint64_t(row[seq & 0xff]) << 16 * 2;
		seq >>= 8;
		b |= uint64_t(row[seq & 0xff]) << 16 * 3;
		seq >>= 8;
		uint64_t c = uint64_t(row[seq & 0xff]);
		seq >>= 8;
		c |= uint64_t(row[seq & 0xff]) << 16;
		seq >>= 8;
		c |= uint64_t(row[seq & 0xff]) << 16 * 2;
		seq >>= 8;
		c |= uint64_t(row[seq & 0xff]) << 16 * 3;
		data_ = _mm_set_epi64x(c, b);
	}

	score_vector(unsigned a, const __m128i &seq, const score_vector &bias)
	{
#ifdef __SSSE3__
		const __m128i *row = reinterpret_cast<const __m128i*>(&score_matrix.matrix8u()[a << 5]);

		__m128i high_mask = _mm_slli_epi16(_mm_and_si128(seq, _mm_set1_epi8('\x10')), 3);
		__m128i seq_low = _mm_or_si128(seq, high_mask);
		__m128i seq_high = _mm_or_si128(seq, _mm_xor_si128(high_mask, _mm_set1_epi8('\x80')));

		__m128i r1 = _mm_load_si128(row);
		__m128i r2 = _mm_load_si128(row + 1);
		__m128i s1 = _mm_shuffle_epi8(r1, seq_low);
		__m128i s2 = _mm_shuffle_epi8(r2, seq_high);
		data_ = _mm_subs_epi16(_mm_and_si128(_mm_or_si128(s1, s2), _mm_set1_epi16(255)), bias.data_);
#endif
	}

	score_vector& zero()
	{
		data_ = _mm_set1_epi16(SHRT_MIN);
		return *this;
	}

	score_vector operator+(const score_vector &rhs) const
	{
		return score_vector(_mm_adds_epi16(data_, rhs.data_));
	}

	score_vector operator-(const score_vector &rhs) const
	{
		return score_vector(_mm_subs_epi16(data_, rhs.data_));
	}

	score_vector& operator-=(const score_vector &rhs)
	{
		data_ = _mm_subs_epi16(data_, rhs.data_);
		return *this;
	}

	int16_t operator [](unsigned i) const
	{
		return *(((int16_t*)&data_) + i); // ^ 0x8000;
	}

	score_vector& max(const score_vector &rhs)
	{
		data_ = _mm_max_epi16(data_, rhs.data_);
		return *this;
	}

	friend score_vector max(const score_vector& lhs, const score_vector &rhs)
	{
		return score_vector(_mm_max_epi16(lhs.data_, rhs.data_));
	}

	uint16_t cmpeq(const score_vector &rhs) const
	{
		return _mm_movemask_epi8(_mm_cmpeq_epi16(data_, rhs.data_));
	}

	__m128i cmpgt(const score_vector &rhs) const
	{
		return _mm_cmpgt_epi16(data_, rhs.data_);
	}

	operator int() const
	{
		return this->operator[](0);
	}

	__m128i data_;

};

#endif

template<typename _t>
struct ScoreTraits
{};

template<>
struct ScoreTraits<int32_t>
{
	enum { CHANNELS = 1 };
	typedef int32_t Score;
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
};

inline void store_sv(int32_t sv, int32_t *dst)
{
	*dst = sv;
}

#ifdef __SSE2__

template<>
struct ScoreTraits<score_vector<int16_t> >
{
	enum { CHANNELS = 8 };
	typedef int16_t Score;
	static score_vector<int16_t> zero()
	{
		return score_vector<int16_t>().zero();
	}
	static void saturate(score_vector<int16_t> &v)
	{
	}
	static int16_t zero_score()
	{
		return SHRT_MIN;
	}
	static int int_score(Score s)
	{
		return (uint16_t)s ^ 0x8000;
	}
	static int16_t max_score()
	{
		return SHRT_MAX;
	}
};

template<typename _t, typename _p>
inline void store_sv(const score_vector<_t> &sv, _p *dst)
{
	_mm_storeu_si128((__m128i*)dst, sv.data_);
}

#endif

#endif /* SCORE_VECTOR_H_ */

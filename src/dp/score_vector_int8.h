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

#ifndef SCORE_VECTOR_INT8_H_
#define SCORE_VECTOR_INT8_H_

#include "score_vector.h"

#ifdef __SSE4_1__

template<>
struct score_vector<int8_t>
{

	score_vector():
#ifdef __AVX2__
		data_(_mm256_set1_epi8(std::numeric_limits<char>::min()))
#else
		data_(_mm_set1_epi8(std::numeric_limits<char>::min()))
#endif
	{}

#ifdef __AVX2__
	explicit score_vector(__m256i data) :
#else
	explicit score_vector(__m128i data):
#endif
		data_(data)
	{}

	explicit score_vector(int8_t x):
#ifdef __AVX2__
		data_(_mm256_set1_epi8(x))
#else
		data_(_mm_set1_epi8(x))
#endif
	{}

	explicit score_vector(int x):
#ifdef __AVX2__
		data_(_mm256_set1_epi8(x))
#else
		data_(_mm_set1_epi8(x))
#endif
	{}

	explicit score_vector(const int8_t* s) :
#ifdef __AVX2__
		data_(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(s)))
#else
		data_(_mm_loadu_si128(reinterpret_cast<const __m128i*>(s)))
#endif
	{ }

	explicit score_vector(const uint8_t* s) :
#ifdef __AVX2__
		data_(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(s)))
#else
		data_(_mm_loadu_si128(reinterpret_cast<const __m128i*>(s)))
#endif
	{ }

#ifdef __AVX2__
	score_vector(unsigned a, __m256i seq)
#else
	score_vector(unsigned a, __m128i seq)
#endif
	{
#ifdef __AVX2__
#elif defined(__SSSE3__)
		const __m128i* row = reinterpret_cast<const __m128i*>(&score_matrix.matrix8()[a << 5]);

		__m128i high_mask = _mm_slli_epi16(_mm_and_si128(seq, _mm_set1_epi8('\x10')), 3);
		__m128i seq_low = _mm_or_si128(seq, high_mask);
		__m128i seq_high = _mm_or_si128(seq, _mm_xor_si128(high_mask, _mm_set1_epi8('\x80')));

		__m128i r1 = _mm_load_si128(row);
		__m128i r2 = _mm_load_si128(row + 1);
		__m128i s1 = _mm_shuffle_epi8(r1, seq_low);
		__m128i s2 = _mm_shuffle_epi8(r2, seq_high);
		data_ = _mm_or_si128(s1, s2);
#else
		const int8_t* row(&score_matrix.matrix8()[a << 5]);
		alignas(16) int8_t data[16], scores[16];
		_mm_store_si128((__mm128i*)data, seq);
		for (unsigned i = 0; i < 16; i++)
			scores[i] = row[data[i]];
		data_ = _mm_load_si128((const __m128i*)scores);
#endif
	}

	score_vector operator+(const score_vector &rhs) const
	{
#ifdef __AVX2__
		return score_vector(_mm256_adds_epi8(data_, rhs.data_));
#else
		return score_vector(_mm_adds_epi8(data_, rhs.data_));
#endif
	}

	score_vector operator-(const score_vector &rhs) const
	{
#ifdef __AVX2__
		return score_vector(_mm256_subs_epi8(data_, rhs.data_));
#else
		return score_vector(_mm_subs_epi8(data_, rhs.data_));
#endif
	}

	score_vector& operator+=(const score_vector& rhs) {
#ifdef __AVX2__
		data_ = _mm256_adds_epi8(data_, rhs.data_);
#else
		data_ = _mm_adds_epi8(data_, rhs.data_);
#endif
		return *this;
	}

	score_vector& operator-=(const score_vector& rhs)
	{
#ifdef __AVX2__
		data_ = _mm256_subs_epi8(data_, rhs.data_);
#else
		data_ = _mm_subs_epi8(data_, rhs.data_);
#endif
		return *this;
	}

	score_vector& operator &=(const score_vector& rhs) {
#ifdef __AVX2__
		data_ = _mm256_and_si256(data_, rhs.data_);
#else
		data_ = _mm_and_si128(data_, rhs.data_);
#endif
		return *this;
	}

	score_vector& operator++() {
#ifdef __AVX2__
		data_ = _mm256_adds_epi8(data_, _mm256_set1_epi8(1));
#else
		data_ = _mm_adds_epi8(data_, _mm_set1_epi8(1));
#endif
		return *this;
	}

#ifdef __AVX2__
	__m256i operator==(const score_vector& rhs) const
#else
	__m128i operator==(const score_vector &rhs) const
#endif
	{
#ifdef __AVX2__
		return _mm256_cmpeq_epi8(data_, rhs.data_);
#else
		return _mm_cmpeq_epi8(data_, rhs.data_);
#endif
	}

	int operator [](unsigned i) const
	{
		return *(((uint8_t*)&data_) + i);
	}

	void set(unsigned i, uint8_t v)
	{
		*(((uint8_t*)&data_) + i) = v;
	}

	score_vector& max(const score_vector &rhs)
	{
#ifdef __AVX2__
		data_ = _mm256_max_epi8(data_, rhs.data_);
#else
		data_ = _mm_max_epi8(data_, rhs.data_);
#endif
		return *this;
	}

	score_vector& min(const score_vector &rhs)
	{
#ifdef __AVX2__
		data_ = _mm256_min_epi8(data_, rhs.data_);
#else
		data_ = _mm_min_epi8(data_, rhs.data_);
#endif
		return *this;
	}

	friend score_vector max(const score_vector& lhs, const score_vector &rhs)
	{
#ifdef __AVX2__
		return score_vector(_mm256_max_epi8(lhs.data_, rhs.data_));
#else
		return score_vector(_mm_max_epi8(lhs.data_, rhs.data_));
#endif
	}

	friend score_vector min(const score_vector& lhs, const score_vector &rhs)
	{
#ifdef __AVX2__
		return score_vector(_mm256_min_epi8(lhs.data_, rhs.data_));
#else
		return score_vector(_mm_min_epi8(lhs.data_, rhs.data_));
#endif
	}

	void store(int8_t *ptr) const
	{
#ifdef __AVX2__
		_mm256_storeu_si256((__m256i*)ptr, data_);
#else
		_mm_storeu_si128((__m128i*)ptr, data_);
#endif
	}

	friend std::ostream& operator<<(std::ostream &s, score_vector v)
	{
		int8_t x[16];
		v.store(x);
		for (unsigned i = 0; i < 16; ++i)
			printf("%3i ", (int)x[i]);
		return s;
	}

#ifdef __AVX2__
	__m256i data_;
#else
	__m128i data_;
#endif

};

template<>
struct ScoreTraits<score_vector<int8_t>>
{
	enum { CHANNELS = 16, BITS = 8 };
	typedef int8_t Score;
	typedef uint8_t Unsigned;
#ifdef __AVX2__
	typedef __m256i Register;
#else
	typedef __m128i Register;
#endif
	static score_vector<int8_t> zero() {
		return score_vector<int8_t>();
	}
	static constexpr int8_t max_score() {
		return SCHAR_MAX;
	}
	static int int_score(int8_t s)
	{
		return (int)s - SCHAR_MIN;
	}
	static constexpr int max_int_score() {
		return SCHAR_MAX - SCHAR_MIN;
	}
	static constexpr int8_t zero_score() {
		return SCHAR_MIN;
	}
	static void saturate(score_vector<int8_t> &v) {}
};

static inline score_vector<int8_t> load_sv(const int8_t *x) {
	return score_vector<int8_t>(x);
}

static inline score_vector<int8_t> load_sv(const uint8_t *x) {
	return score_vector<int8_t>(x);
}

#endif

#endif
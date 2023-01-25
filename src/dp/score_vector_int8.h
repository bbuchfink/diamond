/****
DIAMOND protein aligner
Copyright (C) 2013-2022 Max Planck Society for the Advancement of Science e.V.
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
#include "score_vector.h"

namespace DISPATCH_ARCH {

#if ARCH_ID == 3

template<int DELTA>
struct ScoreVector<int8_t, DELTA>
{

	ScoreVector() :
		data_(_mm512_set1_epi8(DELTA))
	{}

	explicit ScoreVector(__m512i data) :
		data_(data)
	{}

	explicit ScoreVector(int8_t x) :
		data_(_mm512_set1_epi8(x))
	{}

	explicit ScoreVector(int x) :
		data_(_mm512_set1_epi8(x))
	{}

	explicit ScoreVector(const int8_t* s) :
		data_(_mm512_loadu_si512(reinterpret_cast<const __m512i*>(s)))
	{ }

	explicit ScoreVector(const uint8_t* s) :
		data_(_mm512_loadu_si512(reinterpret_cast<const __m512i*>(s)))
	{ }

	ScoreVector(unsigned a, __m512i seq)
	{
		/*const __m256i* row_lo = reinterpret_cast<const __m256i*>(&score_matrix.matrix8_low()[a << 5]);
		const __m256i* row_hi = reinterpret_cast<const __m256i*>(&score_matrix.matrix8_high()[a << 5]);

		__m256i high_mask = _mm256_slli_epi16(_mm256_and_si256(seq, _mm256_set1_epi8('\x10')), 3);
		__m256i seq_low = _mm256_or_si256(seq, high_mask);
		__m256i seq_high = _mm256_or_si256(seq, _mm256_xor_si256(high_mask, _mm256_set1_epi8('\x80')));

		__m256i r1 = _mm512_load_si512(row_lo);
		__m256i r2 = _mm512_load_si512(row_hi);

		__m256i s1 = _mm256_shuffle_epi8(r1, seq_low);
		__m256i s2 = _mm256_shuffle_epi8(r2, seq_high);
		data_ = _mm256_or_si256(s1, s2);*/
	}

	ScoreVector operator+(const ScoreVector& rhs) const
	{
		return ScoreVector(_mm512_adds_epi8(data_, rhs.data_));
	}

	ScoreVector operator-(const ScoreVector& rhs) const
	{
		return ScoreVector(_mm512_subs_epi8(data_, rhs.data_));
	}

	ScoreVector& operator+=(const ScoreVector& rhs) {
		data_ = _mm512_adds_epi8(data_, rhs.data_);
		return *this;
	}

	ScoreVector& operator-=(const ScoreVector& rhs)
	{
		data_ = _mm512_subs_epi8(data_, rhs.data_);
		return *this;
	}

	ScoreVector& operator &=(const ScoreVector& rhs) {
		data_ = _mm512_and_si512(data_, rhs.data_);
		return *this;
	}

	ScoreVector& operator++() {
		data_ = _mm512_adds_epi8(data_, _mm512_set1_epi8(1));
		return *this;
	}

	friend ScoreVector blend(const ScoreVector&v, const ScoreVector&w, const ScoreVector&mask) {
		return ScoreVector(); // (_mm256_blendv_epi8(v.data_, w.data_, mask.data_));
	}

	ScoreVector operator==(const ScoreVector&v) const {
		return ScoreVector(); // ScoreVector(_mm256_cmpeq_epi8(data_, v.data_));
	}

	friend uint32_t cmp_mask(const ScoreVector&v, const ScoreVector&w) {
		return 0; // (uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(v.data_, w.data_));
	}

	int operator [](unsigned i) const
	{
		return *(((uint8_t*)&data_) + i);
	}

	ScoreVector& set(unsigned i, uint8_t v)
	{
		*(((uint8_t*)&data_) + i) = v;
		return *this;
	}

	ScoreVector& max(const ScoreVector& rhs)
	{
		data_ = _mm512_max_epi8(data_, rhs.data_);
		return *this;
	}

	ScoreVector& min(const ScoreVector& rhs)
	{
		data_ = _mm512_min_epi8(data_, rhs.data_);
		return *this;
	}

	friend ScoreVector max(const ScoreVector& lhs, const ScoreVector& rhs)
	{
		return ScoreVector(_mm512_max_epi8(lhs.data_, rhs.data_));
	}

	friend ScoreVector min(const ScoreVector& lhs, const ScoreVector& rhs)
	{
		return ScoreVector(_mm512_min_epi8(lhs.data_, rhs.data_));
	}

	void store(int8_t* ptr) const
	{
		_mm512_storeu_si512((__m512i*)ptr, data_);
	}

	friend std::ostream& operator<<(std::ostream& s, ScoreVector v)
	{
		int8_t x[32];
		v.store(x);
		for (unsigned i = 0; i < 32; ++i)
			printf("%3i ", (int)x[i]);
		return s;
	}

	void expand_from_8bit() {}

	__m512i data_;

};

template<int DELTA>
struct ScoreTraits<ScoreVector<int8_t, DELTA>>
{
	enum { CHANNELS = 64 };
	typedef ::DISPATCH_ARCH::SIMD::Vector<int8_t> Vector;
	typedef int8_t Score;
	typedef uint8_t Unsigned;
	typedef uint32_t Mask;
	struct TraceMask {
		static uint64_t make(uint32_t vmask, uint32_t hmask) {
			return (uint64_t)vmask << 32 | (uint64_t)hmask;
		}
		static uint64_t vmask(int channel) {
			return (uint64_t)1 << (channel + 32);
		}
		static uint64_t hmask(int channel) {
			return (uint64_t)1 << channel;
		}
		uint64_t gap;
		uint64_t open;
	};
	static ScoreVector<int8_t, DELTA> zero() {
		return ScoreVector<int8_t, DELTA>();
	}
	static constexpr int8_t max_score() {
		return SCHAR_MAX;
	}
	static int int_score(int8_t s)
	{
		return (int)s - DELTA;
	}
	static constexpr int max_int_score() {
		return SCHAR_MAX - DELTA;
	}
	static constexpr int8_t zero_score() {
		return DELTA;
	}
	static void saturate(ScoreVector<int8_t, DELTA>& v) {}
};

#elif ARCH_ID == 2

template<int DELTA>
struct ScoreVector<int8_t, DELTA>
{

	ScoreVector() :
		data_(_mm256_set1_epi8(DELTA))
	{}

	explicit ScoreVector(__m256i data) :
		data_(data)
	{}

	explicit ScoreVector(int8_t x) :
		data_(_mm256_set1_epi8(x))
	{}

	explicit ScoreVector(int x) :
		data_(_mm256_set1_epi8(x))
	{}

	explicit ScoreVector(const int8_t* s) :
		data_(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(s)))
	{ }

	explicit ScoreVector(const uint8_t* s) :
		data_(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(s)))
	{ }

	ScoreVector(unsigned a, __m256i seq)
	{
		const __m256i* row_lo = reinterpret_cast<const __m256i*>(&score_matrix.matrix8_low()[a << 5]);
		const __m256i* row_hi = reinterpret_cast<const __m256i*>(&score_matrix.matrix8_high()[a << 5]);

		seq = letter_mask(seq);

		__m256i high_mask = _mm256_slli_epi16(_mm256_and_si256(seq, _mm256_set1_epi8('\x10')), 3);
		__m256i seq_low = _mm256_or_si256(seq, high_mask);
		__m256i seq_high = _mm256_or_si256(seq, _mm256_xor_si256(high_mask, _mm256_set1_epi8('\x80')));

		__m256i r1 = _mm256_load_si256(row_lo);
		__m256i r2 = _mm256_load_si256(row_hi);

		__m256i s1 = _mm256_shuffle_epi8(r1, seq_low);
		__m256i s2 = _mm256_shuffle_epi8(r2, seq_high);
		data_ = _mm256_or_si256(s1, s2);
	}

	ScoreVector operator+(const ScoreVector& rhs) const
	{
		return ScoreVector(_mm256_adds_epi8(data_, rhs.data_));
	}

	ScoreVector operator-(const ScoreVector& rhs) const
	{
		return ScoreVector(_mm256_subs_epi8(data_, rhs.data_));
	}

	ScoreVector& operator+=(const ScoreVector& rhs) {
		data_ = _mm256_adds_epi8(data_, rhs.data_);
		return *this;
	}

	ScoreVector& operator-=(const ScoreVector& rhs)
	{
		data_ = _mm256_subs_epi8(data_, rhs.data_);
		return *this;
	}

	ScoreVector& operator &=(const ScoreVector& rhs) {
		data_ = _mm256_and_si256(data_, rhs.data_);
		return *this;
	}

	ScoreVector& operator++() {
		data_ = _mm256_adds_epi8(data_, _mm256_set1_epi8(1));
		return *this;
	}

	friend ScoreVector blend(const ScoreVector&v, const ScoreVector&w, const ScoreVector&mask) {
		return ScoreVector(_mm256_blendv_epi8(v.data_, w.data_, mask.data_));
	}

	ScoreVector operator==(const ScoreVector&v) const {
		return ScoreVector(_mm256_cmpeq_epi8(data_, v.data_));
	}

	ScoreVector operator>(const ScoreVector& v) const {
		return ScoreVector(_mm256_cmpgt_epi8(data_, v.data_));
	}

	friend uint32_t cmp_mask(const ScoreVector&v, const ScoreVector&w) {
		return (uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(v.data_, w.data_));
	}

	int operator [](unsigned i) const
	{
		return *(((int8_t*)&data_) + i);
	}

	ScoreVector& set(unsigned i, int8_t v)
	{
		//*(((uint8_t*)&data_) + i) = v;
		//data_ = _mm256_insert_epi8(data_, v, i);
		alignas(32) std::array<int8_t, 32> s;
		_mm256_store_si256((__m256i *)s.data(), data_);
		s[i] = v;
		data_ = _mm256_load_si256((__m256i*)s.data());
		return *this;
	}

	ScoreVector& max(const ScoreVector& rhs)
	{
		data_ = _mm256_max_epi8(data_, rhs.data_);
		return *this;
	}

	ScoreVector& min(const ScoreVector& rhs)
	{
		data_ = _mm256_min_epi8(data_, rhs.data_);
		return *this;
	}

	friend ScoreVector max(const ScoreVector& lhs, const ScoreVector& rhs)
	{
		return ScoreVector(_mm256_max_epi8(lhs.data_, rhs.data_));
	}

	friend ScoreVector min(const ScoreVector& lhs, const ScoreVector& rhs)
	{
		return ScoreVector(_mm256_min_epi8(lhs.data_, rhs.data_));
	}

	void store(int8_t* ptr) const
	{
		_mm256_storeu_si256((__m256i*)ptr, data_);
	}

	void store_aligned(int8_t* ptr) const
	{
		_mm256_store_si256((__m256i*)ptr, data_);
	}

	friend std::ostream& operator<<(std::ostream& s, ScoreVector v)
	{
		int8_t x[32];
		v.store(x);
		for (unsigned i = 0; i < 32; ++i)
			printf("%3i ", (int)x[i]);
		return s;
	}

	static ScoreVector load_aligned(const int8_t* x) {
		return ScoreVector(_mm256_load_si256((const __m256i*)x));
	}

	void expand_from_8bit() {}

	__m256i data_;

};

template<int i, int DELTA>
static inline int8_t extract(ScoreVector<int8_t, DELTA> sv) {
	return (int8_t)_mm256_extract_epi8(sv.data_, i);
}

template<int DELTA>
static inline void store_expanded(ScoreVector<int8_t, DELTA> sv, int16_t* dst) {
	const __m256i z = _mm256_setzero_si256();
	const __m256i a = _mm256_permute4x64_epi64(sv.data_, 216);
	__m256i b = _mm256_unpacklo_epi8(a, z);
	__m256i c = _mm256_slli_si256(_mm256_cmpgt_epi8(z, b), 1);
	_mm256_store_si256((__m256i*)dst, _mm256_or_si256(b, c));

	b = _mm256_unpackhi_epi8(a, z);
	c = _mm256_slli_si256(_mm256_cmpgt_epi8(z, b), 1);
	_mm256_store_si256((__m256i*)(dst + 16), _mm256_or_si256(b, c));
}

template<int DELTA>
static inline void store_expanded(ScoreVector<int8_t, DELTA> sv, int8_t* dst) {
	_mm256_store_si256((__m256i*)dst, sv.data_);
}

template<int DELTA>
struct ScoreTraits<ScoreVector<int8_t, DELTA>>
{
	enum { CHANNELS = 32 };
	typedef ::DISPATCH_ARCH::SIMD::Vector<int8_t> Vector;
	typedef int8_t Score;
	typedef uint8_t Unsigned;
	typedef uint32_t Mask;
	struct TraceMask {
		static uint64_t make(uint32_t vmask, uint32_t hmask) {
			return (uint64_t)vmask << 32 | (uint64_t)hmask;
		}
		static uint64_t vmask(int channel) {
			return (uint64_t)1 << (channel + 32);
		}
		static uint64_t hmask(int channel) {
			return (uint64_t)1 << channel;
		}
		uint64_t gap;
		uint64_t open;
	};
	static ScoreVector<int8_t, DELTA> zero() {
		return ScoreVector<int8_t, DELTA>();
	}
	static constexpr int8_t max_score() {
		return SCHAR_MAX;
	}
	static int int_score(int8_t s)
	{
		return (int)s - DELTA;
	}
	static constexpr int max_int_score() {
		return SCHAR_MAX - DELTA;
	}
	static constexpr int8_t zero_score() {
		return DELTA;
	}
	static void saturate(ScoreVector<int8_t, DELTA>& v) {}
};

#elif defined(__SSE4_1__)

template<int DELTA>
struct ScoreVector<int8_t, DELTA>
{

	ScoreVector():
		data_(_mm_set1_epi8(DELTA))
	{}

	explicit ScoreVector(__m128i data):
		data_(data)
	{}

	explicit ScoreVector(int8_t x):
		data_(_mm_set1_epi8(x))
	{}

	explicit ScoreVector(int x):
		data_(_mm_set1_epi8(x))
	{}

	explicit ScoreVector(const int8_t* s) :
		data_(_mm_loadu_si128(reinterpret_cast<const __m128i*>(s)))
	{ }

	explicit ScoreVector(const uint8_t* s) :
		data_(_mm_loadu_si128(reinterpret_cast<const __m128i*>(s)))
	{ }

#ifdef __SSSE3__
	ScoreVector(unsigned a, __m128i seq)
	{
		const __m128i* row = reinterpret_cast<const __m128i*>(&score_matrix.matrix8()[a << 5]);

		seq = letter_mask(seq);

		__m128i high_mask = _mm_slli_epi16(_mm_and_si128(seq, _mm_set1_epi8('\x10')), 3);
		__m128i seq_low = _mm_or_si128(seq, high_mask);
		__m128i seq_high = _mm_or_si128(seq, _mm_xor_si128(high_mask, _mm_set1_epi8('\x80')));

		__m128i r1 = _mm_load_si128(row);
		__m128i r2 = _mm_load_si128(row + 1);
		__m128i s1 = _mm_shuffle_epi8(r1, seq_low);
		__m128i s2 = _mm_shuffle_epi8(r2, seq_high);
		data_ = _mm_or_si128(s1, s2);
	}
#endif

	ScoreVector operator+(const ScoreVector&rhs) const
	{
		return ScoreVector(_mm_adds_epi8(data_, rhs.data_));
	}

	ScoreVector operator-(const ScoreVector&rhs) const
	{
		return ScoreVector(_mm_subs_epi8(data_, rhs.data_));
	}

	ScoreVector& operator+=(const ScoreVector& rhs) {
		data_ = _mm_adds_epi8(data_, rhs.data_);
		return *this;
	}

	ScoreVector& operator-=(const ScoreVector& rhs)
	{
		data_ = _mm_subs_epi8(data_, rhs.data_);
		return *this;
	}

	ScoreVector& operator &=(const ScoreVector& rhs) {
		data_ = _mm_and_si128(data_, rhs.data_);
		return *this;
	}

	ScoreVector& operator++() {
		data_ = _mm_adds_epi8(data_, _mm_set1_epi8(1));
		return *this;
	}

	ScoreVector shift_left(int bytes) {
		return _mm_slli_si128(data_, bytes);
	}

	friend ScoreVector blend(const ScoreVector&v, const ScoreVector&w, const ScoreVector&mask) {
		return ScoreVector(_mm_blendv_epi8(v.data_, w.data_, mask.data_));
	}

	ScoreVector operator==(const ScoreVector&v) const {
		return ScoreVector(_mm_cmpeq_epi8(data_, v.data_));
	}

	ScoreVector operator>(const ScoreVector& v) const {
		return ScoreVector(_mm_cmpgt_epi8(data_, v.data_));
	}

	friend uint32_t cmp_mask(const ScoreVector&v, const ScoreVector&w) {
		return _mm_movemask_epi8(_mm_cmpeq_epi8(v.data_, w.data_));
	}

	int operator [](unsigned i) const
	{
		return *(((uint8_t*)&data_) + i);
	}

	ScoreVector& set(unsigned i, uint8_t v)
	{
		*(((uint8_t*)&data_) + i) = v;
		return *this;
	}

	ScoreVector& max(const ScoreVector&rhs)
	{
		data_ = _mm_max_epi8(data_, rhs.data_);
		return *this;
	}

	ScoreVector& min(const ScoreVector&rhs)
	{
		data_ = _mm_min_epi8(data_, rhs.data_);
		return *this;
	}

	friend ScoreVector max(const ScoreVector& lhs, const ScoreVector&rhs)
	{
		return ScoreVector(_mm_max_epi8(lhs.data_, rhs.data_));
	}

	friend ScoreVector min(const ScoreVector& lhs, const ScoreVector&rhs)
	{
		return ScoreVector(_mm_min_epi8(lhs.data_, rhs.data_));
	}

	void store(int8_t *ptr) const
	{
		_mm_storeu_si128((__m128i*)ptr, data_);
	}

	void store_aligned(int8_t* ptr) const
	{
		_mm_store_si128((__m128i*)ptr, data_);
	}

	friend std::ostream& operator<<(std::ostream &s, ScoreVector v)
	{
		int8_t x[16];
		v.store(x);
		for (unsigned i = 0; i < 16; ++i)
			printf("%3i ", (int)x[i]);
		return s;
	}

	static ScoreVector load_aligned(const int8_t* x) {
		return ScoreVector(_mm_load_si128((const __m128i*)x));
	}

	void expand_from_8bit() {}

	__m128i data_;

};

template<int i, int DELTA>
static inline int8_t extract(ScoreVector<int8_t, DELTA> sv) {
	return 0;
}

template<int DELTA>
struct ScoreTraits<ScoreVector<int8_t, DELTA>>
{
	enum { CHANNELS = 16 };
	typedef ::DISPATCH_ARCH::SIMD::Vector<int8_t> Vector;
	typedef int8_t Score;
	typedef uint8_t Unsigned;
	typedef uint16_t Mask;
	struct TraceMask {
		static uint32_t make(uint32_t vmask, uint32_t hmask) {
			return vmask << 16 | hmask;
		}
		static uint32_t vmask(int channel) {
			return 1 << (channel + 16);
		}
		static uint32_t hmask(int channel) {
			return 1 << channel;
		}
		uint32_t gap;
		uint32_t open;
	};
	static ScoreVector<int8_t, DELTA> zero() {
		return ScoreVector<int8_t, DELTA>();
	}
	static constexpr int8_t max_score() {
		return SCHAR_MAX;
	}
	static int int_score(int8_t s)
	{
		return (int)s - DELTA;
	}
	static constexpr int max_int_score() {
		return SCHAR_MAX - DELTA;
	}
	static constexpr int8_t zero_score() {
		return DELTA;
	}
	static void saturate(ScoreVector<int8_t, DELTA>& v) {}
};

#endif

}

#ifdef __SSE4_1__

template<int DELTA>
static inline int8_t extract_channel(const DISPATCH_ARCH::ScoreVector<int8_t, DELTA>& v, int i) {
	return v[i];
}

template<int DELTA>
static inline void set_channel(DISPATCH_ARCH::ScoreVector<int8_t, DELTA>& v, const int i, const int8_t x) {
	v.set(i, x);
}


#endif
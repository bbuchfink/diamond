/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once
#include <limits.h>
#include "score_vector.h"
#include "stats/score_matrix.h"

namespace DISPATCH_ARCH {

#if ARCH_AVX2_KERNELS

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

	/* One row of a 32x32 score table, prepared for looking up scores by letter.
	   pshufb only indexes 16 entries per 128 bit lane, so the row is split in
	   halves and each half is repeated across both lanes. */
	struct Table {
		explicit Table(const int8_t* row) :
			low(_mm256_broadcastsi128_si256(_mm_loadu_si128(reinterpret_cast<const __m128i*>(row)))),
			high(_mm256_broadcastsi128_si256(_mm_loadu_si128(reinterpret_cast<const __m128i*>(row + 16))))
		{}
		__m256i low, high;
	};

	// Channel k is set to row[letter_mask(seq[k])].
	ScoreVector(const Table& row, const Letter* seq)
	{
		const __m256i s = letter_mask(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(seq)));
		const __m256i high_mask = _mm256_slli_epi16(_mm256_and_si256(s, _mm256_set1_epi8('\x10')), 3);
		const __m256i s1 = _mm256_shuffle_epi8(row.low, _mm256_or_si256(s, high_mask));
		const __m256i s2 = _mm256_shuffle_epi8(row.high, _mm256_or_si256(s, _mm256_xor_si256(high_mask, _mm256_set1_epi8('\x80'))));
		data_ = _mm256_or_si256(s1, s2);
	}

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

// Stores the channels sign extended to 16 bit (unaligned, 32 values).
template<int DELTA>
static inline void store_expanded(ScoreVector<int8_t, DELTA> sv, int16_t* dst) {
	_mm256_storeu_si256((__m256i*)dst, _mm256_cvtepi8_epi16(_mm256_castsi256_si128(sv.data_)));
	_mm256_storeu_si256((__m256i*)(dst + 16), _mm256_cvtepi8_epi16(_mm256_extracti128_si256(sv.data_, 1)));
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

#elif defined(__ARM_NEON)

template<int DELTA>
struct ScoreVector<int8_t, DELTA>
{

	ScoreVector():
		data_(vdupq_n_s8(DELTA))
	{}

	explicit ScoreVector(int8x16_t data):
		data_(data)
	{}

	explicit ScoreVector(int8_t x):
		data_(vdupq_n_s8(x))
	{}

	explicit ScoreVector(int x):
		data_(vdupq_n_s8(x))
	{}

	explicit ScoreVector(const int8_t* s) :
		data_(vld1q_s8(s))
	{ }

	explicit ScoreVector(const uint8_t* s) :
		data_(vreinterpretq_s8_u8(vld1q_u8(s)))
	{ }

#ifdef __aarch64__
	// One row of a 32x32 score table, prepared for looking up scores by letter.
	struct Table {
		explicit Table(const int8_t* row) {
			t.val[0] = vld1q_s8(row);
			t.val[1] = vld1q_s8(row + 16);
		}
		int8x16x2_t t;
	};

	// Channel k is set to row[letter_mask(seq[k])].
	ScoreVector(const Table& row, const Letter* seq) :
		data_(vqtbl2q_s8(row.t, vreinterpretq_u8_s8(letter_mask(vld1q_s8(seq)))))
	{}
#else
	struct Table {
		explicit Table(const int8_t* row) {
			for (int i = 0; i < 4; ++i)
				t.val[i] = vld1_s8(row + 8 * i);
		}
		int8x8x4_t t;
	};

	ScoreVector(const Table& row, const Letter* seq)
	{
		const int8x16_t s = letter_mask(vld1q_s8(seq));
		data_ = vcombine_s8(vtbl4_s8(row.t, vget_low_s8(s)), vtbl4_s8(row.t, vget_high_s8(s)));
	}
#endif

#ifdef __aarch64__
	ScoreVector(unsigned a, int8x16_t seq)
	{
		const int8x16_t* row = reinterpret_cast<const int8x16_t*>(&score_matrix.matrix8()[a << 5]);

		seq = letter_mask(seq);

		int8x16_t high_mask = vreinterpretq_s8_s16(vshlq_n_s16(vreinterpretq_s16_s8(vandq_s8(seq, vdupq_n_s8('\x10'))), 3));
		int8x16_t seq_low   = vorrq_s8(seq, high_mask);
		int8x16_t seq_high  = vorrq_s8(seq, veorq_s8(high_mask, vdupq_n_s8('\x80')));

		int8x16_t r1 = vld1q_s8(reinterpret_cast<const int8_t*>(row));
		int8x16_t r2 = vld1q_s8(reinterpret_cast<const int8_t*>(row + 1));

		int8x16_t s1 = vqtbl1q_s8(r1, vandq_u8(vreinterpretq_u8_s8(seq_low),  vdupq_n_u8(0x8F)));
		int8x16_t s2 = vqtbl1q_s8(r2, vandq_u8(vreinterpretq_u8_s8(seq_high), vdupq_n_u8(0x8F)));
		data_ = vorrq_s8(s1, s2);
	}
#endif

	ScoreVector operator+(const ScoreVector&rhs) const
	{
		return ScoreVector(vqaddq_s8(data_, rhs.data_));
	}

	ScoreVector operator-(const ScoreVector&rhs) const
	{
		return ScoreVector(vqsubq_s8(data_, rhs.data_));
	}

	ScoreVector& operator+=(const ScoreVector& rhs) {
		data_ = vqaddq_s8(data_, rhs.data_);
		return *this;
	}

	ScoreVector& operator-=(const ScoreVector& rhs)
	{
		data_ = vqsubq_s8(data_, rhs.data_);
		return *this;
	}

	ScoreVector& operator &=(const ScoreVector& rhs) {
		data_ = vandq_s8(data_, rhs.data_);
		return *this;
	}

	ScoreVector& operator++() {
		data_ = vqaddq_s8(data_, vdupq_n_s8(1));
		return *this;
	}

	friend ScoreVector blend(const ScoreVector&v, const ScoreVector&w, const ScoreVector&mask) {
		/* Use a signed shift right to create a mask with the sign bit */
		uint8x16_t mask_ = vreinterpretq_u8_s8(vshrq_n_s8(mask.data_, 7));
		return ScoreVector(vbslq_s8(mask_, w.data_, v.data_));
	}

	ScoreVector operator==(const ScoreVector&v) const {
		return ScoreVector(vreinterpretq_s8_u8(vceqq_s8(data_, v.data_)));
	}

	friend uint32_t cmp_mask(const ScoreVector&v, const ScoreVector&w) {
		return vmaskq_s8(vreinterpretq_s8_u8(vceqq_s8(v.data_, w.data_)));
	}

	int operator [](unsigned i) const
	{
		int8_t x[16];
		store(x);
		return x[i];
	}

	ScoreVector& set(const unsigned i, uint8_t v)
	{
		int8_t x[16];
                store(x);
		x[i] = v;
		data_ = vld1q_s8(x);
		return *this;
	}

	ScoreVector& max(const ScoreVector&rhs)
	{
		data_ = vmaxq_s8(data_, rhs.data_);
		return *this;
	}

	ScoreVector& min(const ScoreVector&rhs)
	{
		data_ = vminq_s8(data_, rhs.data_);
		return *this;
	}

	friend ScoreVector max(const ScoreVector& lhs, const ScoreVector&rhs)
	{
		return ScoreVector(vmaxq_s8(lhs.data_, rhs.data_));
	}

	friend ScoreVector min(const ScoreVector& lhs, const ScoreVector&rhs)
	{
		return ScoreVector(vminq_s8(lhs.data_, rhs.data_));
	}

	void store(int8_t *ptr) const
	{
			vst1q_s8(ptr, data_);
	}

	friend std::ostream& operator<<(std::ostream &s, ScoreVector v)
	{
		int8_t x[16];
		v.store(x);
		for (unsigned i = 0; i < 16; ++i)
			printf("%3i ", (int)x[i]);
		return s;
	}

	void expand_from_8bit() {}

	int8x16_t data_;

};

// Stores the channels sign extended to 16 bit (16 values).
template<int DELTA>
static inline void store_expanded(ScoreVector<int8_t, DELTA> sv, int16_t* dst) {
	vst1q_s16(dst, vmovl_s8(vget_low_s8(sv.data_)));
	vst1q_s16(dst + 8, vmovl_s8(vget_high_s8(sv.data_)));
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

	// One row of a 32x32 score table, prepared for looking up scores by letter.
	struct Table {
		explicit Table(const int8_t* row) :
			low(_mm_loadu_si128(reinterpret_cast<const __m128i*>(row))),
			high(_mm_loadu_si128(reinterpret_cast<const __m128i*>(row + 16)))
		{}
		__m128i low, high;
	};

	// Channel k is set to row[letter_mask(seq[k])].
	ScoreVector(const Table& row, const Letter* seq)
	{
		const __m128i s = letter_mask(_mm_loadu_si128(reinterpret_cast<const __m128i*>(seq)));
		const __m128i high_mask = _mm_slli_epi16(_mm_and_si128(s, _mm_set1_epi8('\x10')), 3);
		const __m128i s1 = _mm_shuffle_epi8(row.low, _mm_or_si128(s, high_mask));
		const __m128i s2 = _mm_shuffle_epi8(row.high, _mm_or_si128(s, _mm_xor_si128(high_mask, _mm_set1_epi8('\x80'))));
		data_ = _mm_or_si128(s1, s2);
	}

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

	friend ScoreVector blend(const ScoreVector& v, const ScoreVector& w, const ScoreVector& mask) {
		return ScoreVector(_mm_blendv_epi8(v.data_, w.data_, mask.data_));
	}

	ScoreVector operator==(const ScoreVector& v) const {
		return ScoreVector(_mm_cmpeq_epi8(data_, v.data_));
	}

	ScoreVector operator>(const ScoreVector& v) const {
		return ScoreVector(_mm_cmpgt_epi8(data_, v.data_));
	}

	friend uint32_t cmp_mask(const ScoreVector& v, const ScoreVector& w) {
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

// Stores the channels sign extended to 16 bit (unaligned, 16 values).
template<int DELTA>
static inline void store_expanded(ScoreVector<int8_t, DELTA> sv, int16_t* dst) {
	_mm_storeu_si128((__m128i*)dst, _mm_cvtepi8_epi16(sv.data_));
	_mm_storeu_si128((__m128i*)(dst + 8), _mm_cvtepi8_epi16(_mm_srli_si128(sv.data_, 8)));
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

#if defined(__SSE4_1__) | defined(__ARM_NEON)

template<int DELTA>
static inline int8_t extract_channel(const DISPATCH_ARCH::ScoreVector<int8_t, DELTA>& v, int i) {
	return v[i];
}

template<int DELTA>
static inline void set_channel(DISPATCH_ARCH::ScoreVector<int8_t, DELTA>& v, const int i, const int8_t x) {
	v.set(i, x);
}

#endif
/****
DIAMOND protein aligner
Copyright (C) 2013-2021 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen
						
Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>
Arm NEON port contributed by Martin Larralde <martin.larralde@embl.de>

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

#if ARCH_ID == 2

template<int DELTA>
struct ScoreVector<int8_t, DELTA>
{

	ScoreVector() :
		data_(::SIMD::_mm256_set1_epi8(DELTA))
	{}

	explicit ScoreVector(__m256i data) :
		data_(data)
	{}

	explicit ScoreVector(int8_t x) :
		data_(::SIMD::_mm256_set1_epi8(x))
	{}

	explicit ScoreVector(int x) :
		data_(::SIMD::_mm256_set1_epi8(x))
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

		__m256i high_mask = _mm256_slli_epi16(_mm256_and_si256(seq, ::SIMD::_mm256_set1_epi8('\x10')), 3);
		__m256i seq_low = _mm256_or_si256(seq, high_mask);
		__m256i seq_high = _mm256_or_si256(seq, _mm256_xor_si256(high_mask, ::SIMD::_mm256_set1_epi8('\x80')));

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
		data_ = _mm256_adds_epi8(data_, ::SIMD::_mm256_set1_epi8(1));
		return *this;
	}

	friend ScoreVector blend(const ScoreVector&v, const ScoreVector&w, const ScoreVector&mask) {
		return ScoreVector(_mm256_blendv_epi8(v.data_, w.data_, mask.data_));
	}

	ScoreVector operator==(const ScoreVector&v) const {
		return ScoreVector(_mm256_cmpeq_epi8(data_, v.data_));
	}

	friend uint32_t cmp_mask(const ScoreVector&v, const ScoreVector&w) {
		return (uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(v.data_, w.data_));
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

	friend std::ostream& operator<<(std::ostream& s, ScoreVector v)
	{
		int8_t x[32];
		v.store(x);
		for (unsigned i = 0; i < 32; ++i)
			printf("%3i ", (int)x[i]);
		return s;
	}

	void expand_from_8bit() {}

	__m256i data_;

};

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

#ifdef __AARCH64__
	ScoreVector(unsigned a, int8x16_t seq)
	{
		const int8x16_t* row = reinterpret_cast<const int8x16_t*>(&score_matrix.matrix8()[a << 5]);

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
		/* https://github.com/simd-everywhere/simde/blob/master/simde/x86/sse2.h#L3755 */
		const uint8x16_t MASK = {
			1 << 0, 1 << 1, 1 << 2, 1 << 3, 1 << 4, 1 << 5, 1 << 6, 1 << 7,
			1 << 0, 1 << 1, 1 << 2, 1 << 3, 1 << 4, 1 << 5, 1 << 6, 1 << 7,
		};
		uint8x16_t extended = vceqq_s8(v.data_, w.data_);
		uint8x16_t masked = vandq_u8(MASK, extended);
#ifdef __AARCH64__
		uint8x8x2_t tmp = vzip_u8(vget_low_u8(masked), vget_high_u8(masked));
		uint16x8_t x = vreinterpretq_u16_u8(vcombine_u8(tmp.val[0], tmp.val[1]));
		return vaddvq_u16(x);
#else
		uint16x8_t spliced    = vpadalq_u8(vdupq_n_u16(0), masked);
		uint16x4_t spliced_lo = vget_low_u16(spliced);
		uint16x4_t spliced_hi = vget_high_u16(spliced);
		uint16x8_t zipped = vcombine_u16(spliced_lo, vshl_n_u16(spliced_hi, 8));
		return ::SIMD::vhsumq_u16(zipped);
#endif
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
		data_(::SIMD::_mm_set1_epi8(DELTA))
	{}

	explicit ScoreVector(__m128i data):
		data_(data)
	{}

	explicit ScoreVector(int8_t x):
		data_(::SIMD::_mm_set1_epi8(x))
	{}

	explicit ScoreVector(int x):
		data_(::SIMD::_mm_set1_epi8(x))
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

		__m128i high_mask = _mm_slli_epi16(_mm_and_si128(seq, ::SIMD::_mm_set1_epi8('\x10')), 3);
		__m128i seq_low = _mm_or_si128(seq, high_mask);
		__m128i seq_high = _mm_or_si128(seq, _mm_xor_si128(high_mask, ::SIMD::_mm_set1_epi8('\x80')));

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
		data_ = _mm_adds_epi8(data_, ::SIMD::_mm_set1_epi8(1));
		return *this;
	}

	friend ScoreVector blend(const ScoreVector&v, const ScoreVector&w, const ScoreVector&mask) {
		return ScoreVector(_mm_blendv_epi8(v.data_, w.data_, mask.data_));
	}

	ScoreVector operator==(const ScoreVector&v) const {
		return ScoreVector(_mm_cmpeq_epi8(data_, v.data_));
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

	friend std::ostream& operator<<(std::ostream &s, ScoreVector v)
	{
		int8_t x[16];
		v.store(x);
		for (unsigned i = 0; i < 16; ++i)
			printf("%3i ", (int)x[i]);
		return s;
	}

	void expand_from_8bit() {}

	__m128i data_;

};

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
static inline DISPATCH_ARCH::ScoreVector<int8_t, DELTA> set_channel(const DISPATCH_ARCH::ScoreVector<int8_t, DELTA>& v, const int i, const int8_t x) {
	return DISPATCH_ARCH::ScoreVector<int8_t, DELTA>(v).set(i, x);
}

#endif

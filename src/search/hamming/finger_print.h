/****
Copyright © 2012-2026 Benjamin J. Buchfink <buchfink@gmail.com>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

	http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
****/
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <array>
#include <cstring>
#include "util/simd.h"
#include "basic/value.h"
#include "../search.h"

using std::array;

namespace DISPATCH_ARCH {

#ifdef __AVX512BW__

struct FingerPrint {
    alignas(64) __m512i v;
    static constexpr __mmask64 K48 = (1ULL << 48) - 1;

    explicit FingerPrint(const std::array<char, 48>& a) noexcept
        : v(_mm512_maskz_loadu_epi8(K48, a.data()))
    {}

    static void load(const Letter* q, std::array<char, 48>* dst) noexcept {
#ifdef SEQ_MASK
        __m512i x = _mm512_maskz_loadu_epi8(K48, static_cast<const void*>(q - 16));
		x = _mm512_and_si512(x, _mm512_set1_epi8(LETTER_MASK));
        _mm512_mask_storeu_epi8(dst->data(), K48, x);
#else
        std::copy(q - 16, q + 32, dst->begin());
#endif
    }

    unsigned match(const FingerPrint& rhs) const noexcept {
        __mmask64 m = _mm512_cmpeq_epi8_mask(v, rhs.v) & K48;
		return popcount64(static_cast<unsigned long long>(m));
    }
};

#elif defined(__AVX2__)

struct FingerPrint {
	alignas(32) __m256i v0;
	alignas(16) __m128i v1;

	explicit FingerPrint(const std::array<char, 48>& a) noexcept :
		v0(_mm256_loadu_si256((const __m256i*)a.data())),
		v1(_mm_load_si128((const __m128i*)(a.data() + 32)))
	{
	}

	static void load(const Letter* q, array<char, 48>* dst) noexcept {
#ifdef SEQ_MASK
		__m256i v0 = letter_mask(_mm256_loadu_si256((const __m256i*)(q - 16)));
		__m128i v1 = letter_mask(_mm_loadu_si128((const __m128i*)(q + 16)));
		_mm256_storeu_si256((__m256i*)dst, v0);
		_mm_store_si128((__m128i*)dst + 2, v1);
#else
		std::copy(q - 16, q + 32, dst->begin());
#endif
	}

	static inline int match_block(__m256i a, __m256i b) noexcept {
		return _mm256_movemask_epi8(_mm256_cmpeq_epi8(a, b));
	}

	static int match_block(__m128i x, __m128i y) noexcept {
		return _mm_movemask_epi8(_mm_cmpeq_epi8(x, y));
	}

	unsigned match(const FingerPrint& rhs) const noexcept {
		int m0 = match_block(v0, rhs.v0);
		int m1 = match_block(v1, rhs.v1);
		//uint64_t combined = m0 | (m1 << 32);
		//return popcount64(combined);
		return popcount32(m0) + popcount32(m1);
	}

};

#elif defined(__SSE2__)

struct FingerPrint
{

	explicit FingerPrint(const std::array<char, 48>& a) noexcept :
		r1(_mm_load_si128((const __m128i*)a.data())),
		r2(_mm_load_si128((const __m128i*)(a.data()+16))),
		r3(_mm_load_si128((const __m128i*)(a.data()+32)))
	{
	}

	static void load(const Letter* q, array<char, 48>* dst) noexcept {
#ifdef SEQ_MASK
		__m128i r1 = letter_mask(_mm_loadu_si128((__m128i const*)(q - 16)));
		__m128i r2 = letter_mask(_mm_loadu_si128((__m128i const*)(q)));
		__m128i r3 = letter_mask(_mm_loadu_si128((__m128i const*)(q + 16)));
		_mm_store_si128((__m128i*)dst, r1);
		_mm_store_si128((__m128i*)dst + 1, r2);
		_mm_store_si128((__m128i*)dst + 2, r3);
#else
		std::copy(q - 16, q + 32, dst->begin());
#endif
	}

	static int match_block(__m128i x, __m128i y) noexcept
	{
		return _mm_movemask_epi8(_mm_cmpeq_epi8(x, y));
	}
	unsigned match(const FingerPrint& rhs) const noexcept
	{
		//return popcount64(match_block(r3, rhs.r3) << 32 | match_block(r1, rhs.r1) << 16 | match_block(r2, rhs.r2));
		return popcount32(match_block(r1, rhs.r1) << 16 | match_block(r2, rhs.r2)) + popcount32(match_block(r3, rhs.r3));
	}
	alignas(16) __m128i r1, r2, r3;
};

#elif defined(__ARM_NEON)

struct FingerPrint
{

	explicit FingerPrint(const std::array<char, 48>& a) noexcept :
		r1(vld1q_s8((const signed char*)a.data())),
		r2(vld1q_s8((const signed char*)a.data() + 16)),
		r3(vld1q_s8((const signed char*)a.data() + 32))
	{}

	static void load(const Letter* q, std::array<char, 48>* dst) noexcept {
#ifdef SEQ_MASK
		int8x16_t r1 = letter_mask(vld1q_s8(q - 16));
		int8x16_t r2 = letter_mask(vld1q_s8(q));
		int8x16_t r3 = letter_mask(vld1q_s8(q + 16));
#else
		int8x16_t r1 = vld1q_s8(q - 16);
		int8x16_t r2 = vld1q_s8(q);
		int8x16_t r3 = vld1q_s8(q + 16);
#endif
		vst1q_s8((signed char*)dst->data(), r1);
		vst1q_s8((signed char*)dst->data() + 16, r2);
		vst1q_s8((signed char*)dst->data() + 32, r3);
	}
	
	unsigned match(const FingerPrint& rhs) const noexcept
	{
		const uint8x16_t ONES = vdupq_n_u8(1);
		uint8x16_t s1 = vandq_u8(vceqq_s8(r1, rhs.r1), ONES);
		uint8x16_t s2 = vandq_u8(vceqq_s8(r2, rhs.r2), ONES);
		uint8x16_t s3 = vandq_u8(vceqq_s8(r3, rhs.r3), ONES);
		uint16x8_t acc = vdupq_n_u16(0);
		acc = vpadalq_u8(acc, s1);
		acc = vpadalq_u8(acc, s2);
		acc = vpadalq_u8(acc, s3);
		return vhsumq_u16(acc);
	}

	alignas(16) int8x16_t r1, r2, r3;

};

#else

struct FingerPrint
{
	FingerPrint()
	{
	}
	FingerPrint(const Letter* q) noexcept
	{
		memcpy(r, q - 16, 48);
#ifdef SEQ_MASK
		for (int i = 0; i < 48; ++i)
			r[i] &= LETTER_MASK;
#endif
	}
	unsigned match(const FingerPrint& rhs) const noexcept
	{
		unsigned n = 0;
		for (unsigned i = 0; i < 48; ++i)
			if (r[i] == rhs.r[i])
				++n;
		return n;
	}
	Letter r[48];
};

#endif

template<typename SeedLoc>
static void load_fps(const SeedLoc* p, size_t n, ::Search::Container& v, const SequenceSet& seqs) noexcept
{
	v.resize(n);
	const SeedLoc* end = p + n;
	array<char, 48>* dst = v.data();
	for (; p < end; ++p) {
		FingerPrint::load(seqs.data(*p), dst++);
	}
}

}
/****
DIAMOND protein aligner
Copyright (C) 2016-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
						
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
#include <string>
#include <ostream>
#include "system.h"

#if defined(_M_AMD64) && defined(_MSC_VER)
#define __SSE__
#define __SSE2__
#endif

#ifdef __SSSE3__
#include <tmmintrin.h>
#endif
#ifdef __SSE2__
#include <emmintrin.h>
#endif
#ifdef __SSE4_1__
#include <smmintrin.h>
#endif
#ifdef __AVX2__
#include <immintrin.h>
#endif
#ifdef _MSC_VER
#include <intrin.h>
#endif
#ifdef __ARM_NEON
#include <arm_neon.h>
#endif

namespace SIMD {

enum class Arch { None, Generic, SSE4_1, AVX2, AVX512, NEON };
enum Flags { SSSE3 = 1, POPCNT = 2, SSE4_1 = 4, AVX2 = 8, AVX512 = 16, NEON = 32 };
Arch arch();

std::string features();

}

namespace DISPATCH_ARCH { namespace SIMD {

template<typename T>
struct Vector {};

}}

#ifdef __APPLE__
#ifdef __SSE2__
#define _mm_set1_epi8(v) _mm_set_epi8(v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v)
#define _mm_set1_epi16(v) _mm_set_epi16(v, v, v, v, v, v, v, v)
#endif
#ifdef __AVX2__
#define _mm256_set1_epi8(v) _mm256_set_epi8(v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v)
#define _mm256_set1_epi16(v) _mm256_set_epi16(v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v)
#endif
#endif

#if defined(__GNUC__) && __GNUC__ < 8 && !defined(__clang__)
#define _mm256_set_m128i(v0, v1) _mm256_insertf128_si256(_mm256_castsi128_si256(v1), (v0), 1)
#endif

#ifdef __AVX2__

inline void print_8(__m256i x, std::ostream& s) {
	alignas(32) int8_t v[32];
	_mm256_store_si256((__m256i*)v, x);
	for (unsigned i = 0; i < 32; ++i)
		s << (int)v[i] << ' ';
}

inline void print_16(__m256i x, std::ostream& s) {
	alignas(32) int16_t v[16];
	_mm256_store_si256((__m256i*)v, x);
	for (unsigned i = 0; i < 16; ++i)
		s << (int)v[i] << ' ';
}

#endif

#ifdef __ARM_NEON

inline uint64_t vhsumq_u64(uint64x2_t x) {
#ifdef __aarch64__
	return vaddvq_u64(x);
#else
	return vgetq_lane_u64(x, 0) + vgetq_lane_u64(x, 1);
#endif
}

inline uint64_t vhsumq_u32(uint32x4_t x) {
	uint64x2_t y = vpaddlq_u32(x);
	return vhsumq_u64(y);
}

inline uint64_t vhsumq_u16(uint16x8_t x) {
	uint32x4_t y = vpaddlq_u16(x);
#ifdef __aarch64__
	return vaddvq_u32(y);
#else
	return vhsumq_u32(y);
#endif
}

inline uint64_t vhsumq_u8(uint8x16_t x) {
	uint16x8_t y = vpaddlq_u8(x);
#ifdef __aarch64__
	return vaddvq_u16(y);
#else
	return vhsumq_u16(y);
#endif
}

inline uint16_t vmaskq_s8(int8x16_t x) {
	/* https://github.com/simd-everywhere/simde/blob/master/simde/x86/sse2.h#L3755 */
	const uint8x16_t MASK = {
		1 << 0, 1 << 1, 1 << 2, 1 << 3, 1 << 4, 1 << 5, 1 << 6, 1 << 7,
		1 << 0, 1 << 1, 1 << 2, 1 << 3, 1 << 4, 1 << 5, 1 << 6, 1 << 7,
	};
	uint8x16_t  extended = vreinterpretq_u8_s8(vshrq_n_s8(x, 7));
	uint8x16_t  masked   = vandq_u8(MASK, extended);
	uint8x8x2_t zipped   = vzip_u8(vget_low_u8(masked), vget_high_u8(masked));
	uint16x8_t  spliced  = vreinterpretq_u16_u8(vcombine_u8(zipped.val[0], zipped.val[1]));
#ifdef __aarch64__
	return vaddvq_u16(spliced);
#else
	return ::SIMD::vhsumq_u64(vreinterpretq_u64_u16(spliced));
#endif
}

#endif

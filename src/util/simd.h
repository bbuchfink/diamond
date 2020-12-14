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

#pragma once
#include <string>
#include <functional>
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

namespace SIMD {

enum class Arch { None, Generic, SSE4_1, AVX2 };
enum Flags { SSSE3 = 1, POPCNT = 2, SSE4_1 = 4, AVX2 = 8 };
Arch arch();

#ifdef __SSE__

#define DECL_DISPATCH(ret, name, param) namespace ARCH_GENERIC { ret name param; }\
namespace ARCH_SSE4_1 { ret name param; }\
namespace ARCH_AVX2 { ret name param; }\
inline std::function<decltype(ARCH_GENERIC::name)> dispatch_target_##name() {\
switch(::SIMD::arch()) {\
case ::SIMD::Arch::SSE4_1: return ARCH_SSE4_1::name;\
case ::SIMD::Arch::AVX2: return ARCH_AVX2::name;\
default: return ARCH_GENERIC::name;\
}}\
const std::function<decltype(ARCH_GENERIC::name)> name = dispatch_target_##name();

#else

#define DECL_DISPATCH(ret, name, param) namespace ARCH_GENERIC { ret name param; }\
inline std::function<decltype(ARCH_GENERIC::name)> dispatch_target_##name() {\
return ARCH_GENERIC::name;\
}\
const std::function<decltype(ARCH_GENERIC::name)> name = dispatch_target_##name();

#endif

std::string features();

}

namespace DISPATCH_ARCH { namespace SIMD {

template<typename _t>
struct Vector {};

}}

namespace SIMD {

#ifdef __SSE2__

static inline __m128i _mm_set1_epi8(char v) {
#ifdef __APPLE__
	return _mm_set_epi8(v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v);
#else
	return ::_mm_set1_epi8(v);
#endif
}

static inline __m128i _mm_set1_epi16(short v) {
#ifdef __APPLE__
	return _mm_set_epi16(v, v, v, v, v, v, v, v);
#else
	return ::_mm_set1_epi16(v);
#endif
}

#endif

#ifdef __AVX2__

static inline __m256i _mm256_set1_epi8(char v) {
#ifdef __APPLE__
	return _mm256_set_epi8(v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v);
#else
	return ::_mm256_set1_epi8(v);
#endif
}

static inline __m256i _mm256_set1_epi16(short v) {
#ifdef __APPLE__
	return _mm256_set_epi16(v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v);
#else
	return ::_mm256_set1_epi16(v);
#endif
}

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

}
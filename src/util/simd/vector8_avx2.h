/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <stdint.h>
#include "../simd.h"

namespace DISPATCH_ARCH { namespace SIMD {

template<>
struct Vector<int8_t> {

	static constexpr size_t LANES = 32;

	Vector()
	{}

	Vector(const signed char* p) :
		v(_mm256_loadu_si256((const __m256i*)p))
	{}

	operator __m256i() const {
		return v;
	}

	__m256i v;

};

template<>
struct Vector<int16_t> {

	static constexpr size_t LANES = 16;

	Vector()
	{}

	Vector(const int16_t* p) :
		v(_mm256_loadu_si256((const __m256i*)p))
	{}

	operator __m256i() const {
		return v;
	}

	__m256i v;

};

template<>
struct Vector<int32_t> {

	static constexpr size_t LANES = 1;

	Vector()
	{}

	Vector(const int32_t* p) :
		v(*p)
	{}

	operator int32_t() const {
		return v;
	}

	int32_t v;

};

template<>
struct Traits<float> {
	static constexpr size_t LANES = 8;
	using Register = __m256;
};

static inline __m256 add(__m256 a, __m256 b) {
	return _mm256_add_ps(a, b);
}

static inline __m256 zero(__m256) {
	return _mm256_setzero_ps();
}

static inline __m256 set(float x, __m256) {
	return _mm256_set1_ps(x);
}

static inline __m256 unaligned_load(const float* p, __m256) {
	return _mm256_loadu_ps(p);
}

static inline __m256 load(const float* p, __m256) {
	return _mm256_load_ps(p);
}

static inline void unaligned_store(__m256 v, float* p) {
	_mm256_storeu_ps(p, v);
}

static inline void store(__m256 v, float* p) {
	_mm256_store_ps(p, v);
}

static inline __m256 mul(__m256 a, __m256 b) {
	return _mm256_mul_ps(a, b);
}

static inline __m256 fmadd(__m256 a, __m256 b, __m256 c) {
#if defined(__FMA__) || _MSC_VER
	return _mm256_fmadd_ps(a, b, c);
#else
	return add(mul(a, b), c);
#endif
}

static inline float hsum(__m256 a) {
	__m128 vlow = _mm256_castps256_ps128(a);
	__m128 vhigh = _mm256_extractf128_ps(a, 1);
	__m128 vsum = _mm_add_ps(vlow, vhigh);
	vsum = _mm_hadd_ps(vsum, vsum);
	vsum = _mm_hadd_ps(vsum, vsum);
	return _mm_cvtss_f32(vsum);
}

}}
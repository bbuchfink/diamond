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

	static constexpr size_t LANES = 16;

	Vector()
	{}

	Vector(const signed char *p):
		v(_mm_loadu_si128((const __m128i*)p))
	{}

	operator __m128i() const {
		return v;
	}

	__m128i v;

};

template<>
struct Vector<int16_t> {

	static constexpr size_t LANES = 8;

	Vector()
	{}

	Vector(const int16_t* p) :
		v(_mm_loadu_si128((const __m128i*)p))
	{}

	operator __m128i() const {
		return v;
	}

	__m128i v;

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
	static constexpr size_t LANES = 4;
	using Register = __m128;
};

static inline __m128 add(__m128 a, __m128 b) {
	return _mm_add_ps(a, b);
}

static inline __m128 zero(__m128) {
	return _mm_setzero_ps();
}

static inline __m128 set(float x, __m128) {
	return _mm_set_ss(x);
}

static inline __m128 unaligned_load(const float* p, __m128) {
	return _mm_loadu_ps(p);
}

static inline __m128 load(const float* p, __m128) {
	return _mm_load_ps(p);
}

static inline void unaligned_store(__m128 v, float* p) {
}

static inline void store(__m128 v, float* p) {
}

static inline __m128 mul(__m128 a, __m128 b) {
	return _mm_mul_ps(a, b);
}

static inline __m128 fmadd(__m128 a, __m128 b, __m128 c) {
	return __m128();
}

static inline float hsum(__m128 a) {
	return 0;
}

}}
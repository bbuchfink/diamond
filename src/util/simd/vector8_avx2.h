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
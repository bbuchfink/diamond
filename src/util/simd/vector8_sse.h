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
	return _mm_set_ps1(x);
}

static inline __m128 unaligned_load(const float* p, __m128) {
	return _mm_loadu_ps(p);
}

static inline __m128 load(const float* p, __m128) {
	return _mm_load_ps(p);
}

static inline void unaligned_store(__m128 v, float* p) {
	_mm_storeu_ps(p, v);
}

static inline void store(__m128 v, float* p) {
	_mm_store_ps(p, v);
}

static inline __m128 mul(__m128 a, __m128 b) {
	return _mm_mul_ps(a, b);
}

static inline __m128 fmadd(__m128 a, __m128 b, __m128 c) {
	return add(mul(a, b), c);
}

static inline float hsum(__m128 v) {
	__m128 shuf = _mm_shuffle_ps(v, v, _MM_SHUFFLE(2, 3, 0, 1));
	__m128 sums = _mm_add_ps(v, shuf);
	shuf = _mm_shuffle_ps(sums, sums, _MM_SHUFFLE(1, 0, 3, 2));
	sums = _mm_add_ss(sums, shuf);
	return _mm_cvtss_f32(sums);
}

}}
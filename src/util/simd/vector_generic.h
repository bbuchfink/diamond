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

	static constexpr size_t LANES = 1;

	Vector()
	{}

	Vector(const signed char* p) :
		v(*p)
	{}

	int8_t v;

};

template<>
struct Vector<int16_t> {

	static constexpr size_t LANES = 1;

	Vector()
	{}

	Vector(const int16_t* p) :
		v(*p)
	{}

	int16_t v;

};

template<>
struct Vector<int32_t> {

	static constexpr size_t LANES = 1;

	Vector()
	{}

	Vector(const int32_t* p) :
		v(*p)
	{}

	void store(float* p) const {		
	}

	operator int32_t() const {
		return v;
	}

	int32_t v;

};

template<>
struct Traits<float> {
	static constexpr size_t LANES = 1;
	using Register = float;
};

template<>
struct Traits<float> {
	static constexpr size_t LANES = 1;
	using Register = float;
};

static inline float zero(float) {
	return 0.0f;
}

static inline float set(float x, float) {
	return x;
}

static inline float load(const float* p, float) {
	return *p;
}

static inline float unaligned_load(const float* p, float) {
	return *p;
}

static inline void store(float v, float* p) {
	*p = v;
}

static inline void unaligned_store(float v, float* p) {
	*p = v;
}

static inline float add(float a, float b) {
	return a + b;
}

static inline float mul(float a, float b) {
	return a * b;
}

static inline float fmadd(float a, float b, float c) {
	return a * b + c;
}

static inline float hsum(float a) {
	return a;
}

}}
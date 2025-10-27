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
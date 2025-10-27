/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>
Arm NEON port contributed by Martin Larralde <martin.larralde@embl.de>

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

	static constexpr size_t CHANNELS = 16;

	Vector()
	{}

	Vector(const signed char *p):
		v(vld1q_s8(p))
	{}

	operator int8x16_t() const {
		return v;
	}

	int8x16_t v;

};

template<>
struct Vector<int16_t> {

	static constexpr size_t CHANNELS = 8;

	Vector()
	{}

	Vector(const int16_t* p) :
		v(vld1q_s16(p))
	{}

	operator int16x8_t() const {
		return v;
	}

	int16x8_t v;

};

template<>
struct Vector<int32_t> {

	static constexpr size_t CHANNELS = 1;

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
	using Register = float32x4_t;
};

static inline float32x4_t add(float32x4_t a, float32x4_t b) {
	return vaddq_f32(a, b);
}

static inline float32x4_t zero(float32x4_t) {
	return vdupq_n_f32(0.0f);
}

// Sets only lane 0 to x (like _mm_set_ss), other lanes = 0
static inline float32x4_t set(float x, float32x4_t) {
	float32x4_t z = vdupq_n_f32(0.0f);
	return vsetq_lane_f32(x, z, 0);
}

static inline float32x4_t unaligned_load(const float* p, float32x4_t) {
	// NEON uses the same intrinsic for aligned/unaligned loads
	return vld1q_f32(p);
}

static inline float32x4_t load(const float* p, float32x4_t) {
	return vld1q_f32(p);
}

static inline void unaligned_store(float32x4_t v, float* p) {
	vst1q_f32(p, v);
}

static inline void store(float32x4_t v, float* p) {
	vst1q_f32(p, v);
}

static inline float32x4_t mul(float32x4_t a, float32x4_t b) {
	return vmulq_f32(a, b);
}

static inline float32x4_t fmadd(float32x4_t a, float32x4_t b, float32x4_t c) {
	// c + a*b  (matches add(mul(a,b), c))
#if defined(__aarch64__)
	return vfmaq_f32(c, a, b);
#else
	// vmlaq_f32 is fused on CPUs with FMA; otherwise it may lower to mul+add
	return vmlaq_f32(c, a, b);
#endif
}

static inline float hsum(float32x4_t v) {
#if defined(__aarch64__)
	return vaddvq_f32(v);
#else
	// pairwise reduce: (v0+v1) + (v2+v3)
	float32x2_t vlow = vget_low_f32(v);
	float32x2_t vhigh = vget_high_f32(v);
	float32x2_t sum2 = vadd_f32(vlow, vhigh);      // [v0+v2, v1+v3]
	float32x2_t sum1 = vpadd_f32(sum2, sum2);      // [ (v0+v2)+(v1+v3), ... ]
	return vget_lane_f32(sum1, 0);
#endif
}

}}
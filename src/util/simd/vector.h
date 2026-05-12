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
#include "../simd.h"

#if ARCH_ID == 3
#include "vector8_avx512.h"
#elif ARCH_ID == 2
#include "vector8_avx2.h"
#elif defined(__ARM_NEON)
#include "vector8_neon.h"
#elif defined(__SSE2__)
#include "vector8_sse.h"
#else
#include "vector_generic.h"
#endif

namespace DISPATCH_ARCH { namespace SIMD {
	
static inline float sum(const float* x, size_t n) {
	using Register = Traits<float>::Register;
	constexpr size_t L = Traits<float>::LANES;
	size_t i = 0;
	Register acc = zero(Register());
	for (; i + L <= n; i += L) {
		acc = add(acc, unaligned_load(x + i, Register()));
	}
	float sum = hsum(acc);
	for (; i < n; ++i) sum += x[i];
	return sum;
}

/*static inline void set(float* dst, float val, size_t n) {
	using Register = Traits<float>::Register;
	constexpr size_t L = Traits<float>::LANES;
	const Register s = set(val, Register());
	size_t i = 0;
	for (; i + L <= n; i += L) unaligned_store(s, dst + i);
	for (; i < n; ++i) dst[i] = val;
}*/

static inline void scale(float* dst, float scale, size_t n) {
	using Register = Traits<float>::Register;
	constexpr size_t L = Traits<float>::LANES;
	const Register s = set(scale, Register());
	size_t i = 0;
	for (; i + L <= n; i += L)
		unaligned_store(mul(unaligned_load(dst + i, Register()), s), dst + i);
	for (; i < n; ++i) dst[i] *= scale;
}

/*static inline void add(float* dst, const float* src, size_t n) {
	using Register = Traits<float>::Register;
	constexpr size_t L = Traits<float>::LANES;
	size_t i = 0;
	for (; i + L <= n; i += L) {
		Register a = unaligned_load(dst + i, Register());
		Register b = unaligned_load(src + i, Register());
		unaligned_store(add(a, b), dst + i);
	}
	for (; i < n; ++i) dst[i] += src[i];
}

static inline void mul(float* dst, const float* src, size_t n) {
	using Register = Traits<float>::Register;
	constexpr size_t L = Traits<float>::LANES;
	size_t i = 0;
	for (; i + L <= n; i += L) {
		Register a = unaligned_load(dst + i, Register());
		Register b = unaligned_load(src + i, Register());
		unaligned_store(mul(a, b), dst + i);
	}
	for (; i < n; ++i) dst[i] *= src[i];
}

static inline void scale_to(float* dst, const float* src, float scale, size_t n) {
	using Register = Traits<float>::Register;
	constexpr size_t L = Traits<float>::LANES;
	const Register s = set(scale, Register());
	size_t i = 0;
	for (; i + L <= n; i += L) {
		const Register a = unaligned_load(src + i, Register());
		unaligned_store(mul(a, s), dst + i);
	}
	for (; i < n; ++i) dst[i] = src[i] * scale;
}

static inline void mul_to(float* dst, const float* a, const float* b, size_t n) {
	using Register = Traits<float>::Register;
	constexpr size_t L = Traits<float>::LANES;
	size_t i = 0;
	for (; i + L <= n; i += L) {
		Register va = unaligned_load(a + i, Register());
		Register vb = unaligned_load(b + i, Register());
		unaligned_store(mul(va, vb), dst + i);
	}
	for (; i < n; ++i) dst[i] = a[i] * b[i];
}

static inline void add_scalar_inplace(float* dst, float s, size_t n) {
	using Register = Traits<float>::Register;
	constexpr size_t L = Traits<float>::LANES;
	const Register v = set(s, Register());
	size_t i = 0;
	for (; i + L <= n; i += L) {
		Register a = unaligned_load(dst + i, Register());
		unaligned_store(add(a, v), dst + i);
	}
	for (; i < n; ++i) dst[i] += s;
}*/

}}
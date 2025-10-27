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
#include "../simd.h"

#if ARCH_ID == 3
#include "vector8_avx512.h"
#elif ARCH_ID == 2
#include "vector8_avx2.h"
#elif ARCH_ID == 4
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
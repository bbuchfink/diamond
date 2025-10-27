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

/* Based on tantan by Martin C. Frith. See:
http://cbrc3.cbrc.jp/~martin/tantan/
A new repeat-masking method enables specific detection of homologous sequences, MC Frith, Nucleic Acids Research 2011 39(4):e23. */

#include <vector>
#include "../basic/value.h"
#include "def.h"
#include "../util/simd/dispatch.h"
#include "util/simd/vector.h"
#include "masking.h"

using std::fill;

namespace Util { namespace tantan { namespace DISPATCH_ARCH {

static inline float forward_step(
    float* __restrict f, const float* __restrict d, const float* __restrict e_seg,
    float& b, float f2f, float p_repeat_end, float b2b, float f_sum_prev)
{
	using namespace ::DISPATCH_ARCH::SIMD;
	using Register = Traits<float>::Register;
	constexpr size_t L = Traits<float>::LANES;
    const float b_old = b;
    const Register vf2f   = set(f2f, Register());
    const Register vb_old = set(b_old, Register());

    float f_sum_new = 0.0f;

    for (int off = 0; off < 48; off += L) {
        Register vf = load(f + off, Register());
        const Register vd = load(d + off, Register());
        const Register ve = unaligned_load(e_seg + off, Register());
		Register tmp = fmadd(vf, vf2f, mul(vb_old, vd));
		vf = mul(tmp, ve);
		store(vf, f + off);
		f_sum_new += hsum(vf);
    }
    for (int off = 48; off < 50; ++off) {
        float vf = f[off];
        vf = (vf * f2f + b_old * d[off]) * e_seg[off];
        f[off] = vf;
        f_sum_new += vf;
    }

    b = b_old * b2b + f_sum_prev * p_repeat_end;
    return f_sum_new;
}

static inline float backward_step(
	float* __restrict f, const float* __restrict d, const float* __restrict e_seg,
	float& b, float f2f, float p_repeat_end, float b2b)
{
	using namespace ::DISPATCH_ARCH::SIMD;
	using Register = Traits<float>::Register;
	constexpr size_t L = Traits<float>::LANES;
	const Register vf2f = set(f2f, Register());
	const Register vC = set(p_repeat_end * b, Register());

	float tsum = 0.0f;

	for (int off = 0; off < 48; off += L) {
		Register vf = load(f + off, Register());
		const Register ve = unaligned_load(e_seg + off, Register());
		const Register vd = load(d + off, Register());
		vf = mul(vf, ve);
		Register vt = mul(vf, vd);
		vf = fmadd(vf, vf2f, vC);
		store(vf, f + off);
		tsum += hsum(vt);
	}

	for (int off = 48; off < 50; ++off) {
		float vf = f[off] * e_seg[off];
		tsum += vf * d[off];
		vf = vf * f2f + p_repeat_end * b;
		f[off] = vf;
	}

	b = b2b * b + tsum;
	return tsum;
}

Mask::Ranges mask(
	Letter *seq,
	int len,
	const float **likelihood_ratio_matrix,
	float p_repeat,
	float p_repeat_end,
	float repeat_growth,
	float p_mask,
	int mask_mode)
{
	constexpr int WINDOW  = 50;
	constexpr int RESERVE = 50000;

	Mask::Ranges ranges;
	if (len == 0) return ranges;

	alignas(32) float f[WINDOW];
	alignas(32) float d[WINDOW];

	const float b2b  = 1.0f - p_repeat;
	const float f2f  = 1.0f - p_repeat_end;
	const float b2f0 = p_repeat * (1.0f - repeat_growth) / (1.0f - std::pow(repeat_growth, (float)WINDOW));

	d[WINDOW - 1] = b2f0;
	for (int i = WINDOW - 2; i >= 0; --i)
		d[i] = d[i + 1] * repeat_growth;

	fill(f, f + WINDOW, 0.0f);
	
	thread_local std::vector<float> pb;
	thread_local std::vector<float> scale;
	thread_local std::vector<float> e[AMINO_ACID_COUNT];

	const int pb_cap = std::max(len, RESERVE);
	const int sc_len = std::max((len - 1) / 16 + 1, (RESERVE - 1) / 16 + 1);

	pb.resize(pb_cap);
	scale.resize(sc_len);

	for (size_t aa = 0; aa < (size_t)AMINO_ACID_COUNT; ++aa) {
		auto &E = e[aa];
		const int need = std::max(RESERVE, len + WINDOW);
		E.resize(need);

		const float* L = likelihood_ratio_matrix[aa];
		float* p = &E[len - 1];
		for (int j = 0; j < len; ++j) {
			const uint8_t idx = static_cast<uint8_t>(letter_mask(seq[j]));
			*(p--) = L[(size_t)idx];
		}
		fill(E.data() + len, E.data() + len + WINDOW, (float)0.0);
	}

	float b = 1.0f;
	float f_sum = 0.0f;

	for (int i = 0; i < len; ++i) {
		const uint8_t ltr = static_cast<uint8_t>(letter_mask(seq[i]));
		const float* e_seg = &e[ltr][len - i];
		float f_sum_new = forward_step(f, d, e_seg, b, f2f, p_repeat_end, b2b, f_sum);
		f_sum = f_sum_new;
		if ((i & 15) == 15) {
			const float s = 1.0f / b;
			scale[(size_t)i / 16] = s;
			b *= s;
			::DISPATCH_ARCH::SIMD::scale(f, s, WINDOW);
			f_sum *= s;
		}
		pb[i] = b;
	}

	const float z = b * b2b + ::DISPATCH_ARCH::SIMD::sum(f, WINDOW) * p_repeat_end;
	const float zinv = 1.0f / z;

	b = b2b;
	fill(f, f + WINDOW, p_repeat_end);
	const Letter mask = value_traits.mask_char;
	
	for (int i = len - 1; i >= 0; --i) {
		const float pf = 1.0f - (pb[i] * b * zinv);

		if ((i & 15) == 15) {
			const float s = scale[(size_t)i / 16];
			b *= s;
			::DISPATCH_ARCH::SIMD::scale(f, s, WINDOW);
		}

		const uint8_t ltr = static_cast<uint8_t>(letter_mask(seq[i]));
		const float* e_seg = &e[ltr][len - i];
		backward_step(f, d, e_seg, b, f2f, p_repeat_end, b2b);

		if (pf >= p_mask) {
			if (mask_mode == 1)
				seq[i] = mask;
			else if (mask_mode == 2)
				seq[i] |= Masking::bit_mask;
			ranges.push_front(i);
		}
	}

	return ranges;
}

}

DISPATCH_8(Mask::Ranges, mask, Letter*, seq, int, len, const float**, likelihood_ratio_matrix, float, p_repeat, float, p_repeat_end, float, repeat_decay, float, p_mask, int, mask_mode)

}}
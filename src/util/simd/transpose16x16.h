/****
Copyright � 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>
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
#include <stddef.h>
#include "../simd.h"

#if defined(__SSE2__)

#define UNPACK128_LO_HI_EPI8(a, b) t = r##a; r##a = _mm_unpacklo_epi8(t, r##b); r##b = _mm_unpackhi_epi8(t, r##b);
#define UNPACK128_LO_HI_EPI16(a, b) t = r##a; r##a = _mm_unpacklo_epi16(t, r##b); r##b = _mm_unpackhi_epi16(t, r##b);
#define UNPACK128_LO_HI_EPI32(a, b) t = r##a; r##a = _mm_unpacklo_epi32(t, r##b); r##b = _mm_unpackhi_epi32(t, r##b);
#define UNPACK128_LO_HI_EPI64(a, b) t = r##a; r##a = _mm_unpacklo_epi64(t, r##b); r##b = _mm_unpackhi_epi64(t, r##b);

#define UNPACK256_LO_HI_EPI16(a, b) t = r##a; r##a = _mm256_unpacklo_epi16(t, r##b); r##b = _mm256_unpackhi_epi16(t, r##b);
#define UNPACK256_LO_HI_EPI32(a, b) t = r##a; r##a = _mm256_unpacklo_epi32(t, r##b); r##b = _mm256_unpackhi_epi32(t, r##b);
#define UNPACK256_LO_HI_EPI64(a, b) t = r##a; r##a = _mm256_unpacklo_epi64(t, r##b); r##b = _mm256_unpackhi_epi64(t, r##b);
#define UNPACK256_LO_HI_EPI128(a, b) t = r##a; r##a = _mm256_permute2x128_si256(t, r##b, 32); r##b = _mm256_permute2x128_si256(t, r##b, 49);

static inline void transpose(const signed char **data, size_t n, signed char *out, const __m128i&) {
	__m128i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15, t;
	r0 = r1 = r2 = r3 = r4 = r5 = r6 = r7 = r8 = r9 = r10 = r11 = r12 = r13 = r14 = r15 = _mm_setzero_si128();
	
	switch (n) {
	case 16:
		r0 = _mm_loadu_si128((const __m128i*)*(data++));
	case 15:
		r1 = _mm_loadu_si128((const __m128i*)*(data++));
	case 14:
		r2 = _mm_loadu_si128((const __m128i*)*(data++));
	case 13:
		r3 = _mm_loadu_si128((const __m128i*)*(data++));
	case 12:
		r4 = _mm_loadu_si128((const __m128i*)*(data++));
	case 11:
		r5 = _mm_loadu_si128((const __m128i*)*(data++));
	case 10:
		r6 = _mm_loadu_si128((const __m128i*)*(data++));
	case 9:
		r7 = _mm_loadu_si128((const __m128i*)*(data++));
	case 8:
		r8 = _mm_loadu_si128((const __m128i*)*(data++));
	case 7:
		r9 = _mm_loadu_si128((const __m128i*)*(data++));
	case 6:
		r10 = _mm_loadu_si128((const __m128i*)*(data++));
	case 5:
		r11 = _mm_loadu_si128((const __m128i*)*(data++));
	case 4:
		r12 = _mm_loadu_si128((const __m128i*)*(data++));
	case 3:
		r13 = _mm_loadu_si128((const __m128i*)*(data++));
	case 2:
		r14 = _mm_loadu_si128((const __m128i*)*(data++));
	case 1:
		r15 = _mm_loadu_si128((const __m128i*)*data);
	}

	UNPACK128_LO_HI_EPI8(0, 1)
	UNPACK128_LO_HI_EPI8(2, 3)
	UNPACK128_LO_HI_EPI8(4, 5)
	UNPACK128_LO_HI_EPI8(6, 7)
	UNPACK128_LO_HI_EPI8(8, 9)
	UNPACK128_LO_HI_EPI8(10, 11)
	UNPACK128_LO_HI_EPI8(12, 13)
	UNPACK128_LO_HI_EPI8(14, 15)

	UNPACK128_LO_HI_EPI16(0, 2)
	UNPACK128_LO_HI_EPI16(4, 6)
	UNPACK128_LO_HI_EPI16(8, 10)
	UNPACK128_LO_HI_EPI16(12, 14)
	UNPACK128_LO_HI_EPI16(1, 3)
	UNPACK128_LO_HI_EPI16(5, 7)
	UNPACK128_LO_HI_EPI16(9, 11)
	UNPACK128_LO_HI_EPI16(13, 15)

	UNPACK128_LO_HI_EPI32(0, 4)
	UNPACK128_LO_HI_EPI32(8, 12)
	UNPACK128_LO_HI_EPI32(2, 6)
	UNPACK128_LO_HI_EPI32(10, 14)
	UNPACK128_LO_HI_EPI32(1, 5)
	UNPACK128_LO_HI_EPI32(9, 13)
	UNPACK128_LO_HI_EPI32(3, 7)
	UNPACK128_LO_HI_EPI32(11, 15)

	UNPACK128_LO_HI_EPI64(0, 8)
	UNPACK128_LO_HI_EPI64(4, 12)
	UNPACK128_LO_HI_EPI64(2, 10)
	UNPACK128_LO_HI_EPI64(6, 14)
	UNPACK128_LO_HI_EPI64(1, 9)
	UNPACK128_LO_HI_EPI64(5, 13)
	UNPACK128_LO_HI_EPI64(3, 11)
	UNPACK128_LO_HI_EPI64(7, 15)

	__m128i* ptr = (__m128i*)out;
	_mm_store_si128(ptr++, r0);
	_mm_store_si128(ptr++, r8);
	_mm_store_si128(ptr++, r4);
	_mm_store_si128(ptr++, r12);
	_mm_store_si128(ptr++, r2);
	_mm_store_si128(ptr++, r10);
	_mm_store_si128(ptr++, r6);
	_mm_store_si128(ptr++, r14);
	_mm_store_si128(ptr++, r1);
	_mm_store_si128(ptr++, r9);
	_mm_store_si128(ptr++, r5);
	_mm_store_si128(ptr++, r13);
	_mm_store_si128(ptr++, r3);
	_mm_store_si128(ptr++, r11);
	_mm_store_si128(ptr++, r7);
	_mm_store_si128(ptr, r15);
}

#elif defined(__ARM_NEON)

#define UNPACK128_LO_HI_EPI8(a, b) t_s8 = vtrnq_s8(r##a, r##b); \
        r##a = t_s8.val[0]; \
        r##b = t_s8.val[1];
#define UNPACK128_LO_HI_EPI16(a, b) t_s16 = vtrnq_s16(vreinterpretq_s16_s8(r##a), vreinterpretq_s16_s8(r##b)); \
        r##a = vreinterpretq_s8_s16(t_s16.val[0]); \
        r##b = vreinterpretq_s8_s16(t_s16.val[1]);
#define UNPACK128_LO_HI_EPI32(a, b) t_s32 = vtrnq_s32(vreinterpretq_s32_s8(r##a), vreinterpretq_s32_s8(r##b)); \
        r##a = vreinterpretq_s8_s32(t_s32.val[0]); \
        r##b = vreinterpretq_s8_s32(t_s32.val[1]);
#define STORE_LOW(dst, a, b) vst1q_s8(dst, vcombine_s8(vget_low_s8(r##a), vget_low_s8(r##b)));
#define STORE_HIGH(dst, a, b) vst1q_s8(dst, vcombine_s8(vget_high_s8(r##a), vget_high_s8(r##b)));

static inline void transpose(const signed char **data, size_t n, signed char *out, const int8x16_t&) {
	int8x16_t r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15;
        int8x16x4_t b0, b1, b2, b3;
        int8x16x2_t t_s8;
        int16x8x2_t t_s16;
        int32x4x2_t t_s32;

	r0 = r1 = r2 = r3 = r4 = r5 = r6 = r7 = r8 = r9 = r10 = r11 = r12 = r13 = r14 = r15 = vdupq_n_s8(0);
	
	switch (n) {
	case 16: r0 = vld1q_s8((const int8_t*)*(data++));
	case 15: r1 = vld1q_s8((const int8_t*)*(data++));
	case 14: r2 = vld1q_s8((const int8_t*)*(data++));
	case 13: r3 = vld1q_s8((const int8_t*)*(data++));
	case 12: r4 = vld1q_s8((const int8_t*)*(data++));
	case 11: r5 = vld1q_s8((const int8_t*)*(data++));
	case 10: r6 = vld1q_s8((const int8_t*)*(data++));
	case 9:  r7 = vld1q_s8((const int8_t*)*(data++));
	case 8:  r8 = vld1q_s8((const int8_t*)*(data++));
	case 7:  r9 = vld1q_s8((const int8_t*)*(data++));
	case 6: r10 = vld1q_s8((const int8_t*)*(data++));
	case 5: r11 = vld1q_s8((const int8_t*)*(data++));
	case 4: r12 = vld1q_s8((const int8_t*)*(data++));
	case 3: r13 = vld1q_s8((const int8_t*)*(data++));
	case 2: r14 = vld1q_s8((const int8_t*)*(data++));
	case 1: r15 = vld1q_s8((const int8_t*)*data);
	}

	UNPACK128_LO_HI_EPI8(0, 1)
	UNPACK128_LO_HI_EPI8(2, 3)
	UNPACK128_LO_HI_EPI8(4, 5)
	UNPACK128_LO_HI_EPI8(6, 7)
	UNPACK128_LO_HI_EPI8(8, 9)
	UNPACK128_LO_HI_EPI8(10, 11)
	UNPACK128_LO_HI_EPI8(12, 13)
	UNPACK128_LO_HI_EPI8(14, 15)

	UNPACK128_LO_HI_EPI16(0, 2)
	UNPACK128_LO_HI_EPI16(4, 6)
	UNPACK128_LO_HI_EPI16(8, 10)
	UNPACK128_LO_HI_EPI16(12, 14)
	UNPACK128_LO_HI_EPI16(1, 3)
	UNPACK128_LO_HI_EPI16(5, 7)
	UNPACK128_LO_HI_EPI16(9, 11)
	UNPACK128_LO_HI_EPI16(13, 15)

	UNPACK128_LO_HI_EPI32(0, 4)
	UNPACK128_LO_HI_EPI32(8, 12)
	UNPACK128_LO_HI_EPI32(2, 6)
	UNPACK128_LO_HI_EPI32(10, 14)
	UNPACK128_LO_HI_EPI32(1, 5)
	UNPACK128_LO_HI_EPI32(9, 13)
	UNPACK128_LO_HI_EPI32(3, 7)
	UNPACK128_LO_HI_EPI32(11, 15)

	STORE_LOW(&out[0x00], 0, 8);
	STORE_LOW(&out[0x10], 1, 9);
	STORE_LOW(&out[0x20], 2, 10);
	STORE_LOW(&out[0x30], 3, 11);
	STORE_LOW(&out[0x40], 4, 12);
	STORE_LOW(&out[0x50], 5, 13);
	STORE_LOW(&out[0x60], 6, 14);
	STORE_LOW(&out[0x70], 7, 15);

	STORE_HIGH(&out[0x80], 0, 8);
	STORE_HIGH(&out[0x90], 1, 9);
	STORE_HIGH(&out[0xa0], 2, 10);
	STORE_HIGH(&out[0xb0], 3, 11);
	STORE_HIGH(&out[0xc0], 4, 12);
	STORE_HIGH(&out[0xd0], 5, 13);
	STORE_HIGH(&out[0xe0], 6, 14);
	STORE_HIGH(&out[0xf0], 7, 15);
}

#endif
#ifdef __AVX2__

static inline void transpose_offset(const int16_t** data, size_t n, ptrdiff_t offset, int16_t* out, __m256i) {
	__m256i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15, t;

	r0 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	r1 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	r2 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	r3 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	r4 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	r5 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	r6 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	r7 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	r8 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	r9 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	r10 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	r11 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	r12 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	r13 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	r14 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	r15 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	
	UNPACK256_LO_HI_EPI16(0, 1)
	UNPACK256_LO_HI_EPI16(2, 3)
	UNPACK256_LO_HI_EPI16(4, 5)
	UNPACK256_LO_HI_EPI16(6, 7)
	UNPACK256_LO_HI_EPI16(8, 9)
	UNPACK256_LO_HI_EPI16(10, 11)
	UNPACK256_LO_HI_EPI16(12, 13)
	UNPACK256_LO_HI_EPI16(14, 15)

	UNPACK256_LO_HI_EPI32(0, 2)
	UNPACK256_LO_HI_EPI32(1, 3)
	UNPACK256_LO_HI_EPI32(4, 6)
	UNPACK256_LO_HI_EPI32(5, 7)
	UNPACK256_LO_HI_EPI32(8, 10)
	UNPACK256_LO_HI_EPI32(9, 11)
	UNPACK256_LO_HI_EPI32(12, 14)
	UNPACK256_LO_HI_EPI32(13, 15)

	UNPACK256_LO_HI_EPI64(0, 4)
	UNPACK256_LO_HI_EPI64(2, 6)
	UNPACK256_LO_HI_EPI64(1, 5)
	UNPACK256_LO_HI_EPI64(3, 7)
	UNPACK256_LO_HI_EPI64(8, 12)
	UNPACK256_LO_HI_EPI64(10, 14)
	UNPACK256_LO_HI_EPI64(9, 13)
	UNPACK256_LO_HI_EPI64(11, 15)

	UNPACK256_LO_HI_EPI128(0, 8)
	UNPACK256_LO_HI_EPI128(4, 12)
	UNPACK256_LO_HI_EPI128(2, 10)
	UNPACK256_LO_HI_EPI128(6, 14)
	UNPACK256_LO_HI_EPI128(1, 9)
	UNPACK256_LO_HI_EPI128(5, 13)
	UNPACK256_LO_HI_EPI128(3, 11)
	UNPACK256_LO_HI_EPI128(7, 15)
	
	__m256i* ptr = (__m256i*)out;
	_mm256_store_si256(ptr++, r0);
	_mm256_store_si256(ptr++, r4);
	_mm256_store_si256(ptr++, r2);
	_mm256_store_si256(ptr++, r6);
	_mm256_store_si256(ptr++, r1);
	_mm256_store_si256(ptr++, r5);
	_mm256_store_si256(ptr++, r3);
	_mm256_store_si256(ptr++, r7);
	_mm256_store_si256(ptr++, r8);
	_mm256_store_si256(ptr++, r12);
	_mm256_store_si256(ptr++, r10);
	_mm256_store_si256(ptr++, r14);
	_mm256_store_si256(ptr++, r9);
	_mm256_store_si256(ptr++, r13);
	_mm256_store_si256(ptr++, r11);
	_mm256_store_si256(ptr++, r15);
}

#define LOAD(r) r = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
#define ALIGN(r) r = _mm256_unpacklo_epi8(_mm256_permute4x64_epi64(r, 2),z);

static inline void transpose_offset_8bit(const int16_t** data, size_t n, ptrdiff_t offset, int16_t* out, __m256i) {
	__m256i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15, t;
	//__m128i s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15;
	__m256i z = _mm256_setzero_si256();

	LOAD(r0)
	LOAD(r1)
	LOAD(r2)
	LOAD(r3)
	LOAD(r4)
	LOAD(r5)
	LOAD(r6)
	LOAD(r7)
	LOAD(r8)
	LOAD(r9)
	LOAD(r10)
	LOAD(r11)
	LOAD(r12)
	LOAD(r13)
	LOAD(r14)
	LOAD(r15)

	ALIGN(r0)
	ALIGN(r1)
	ALIGN(r2)
	ALIGN(r3)
	ALIGN(r4)
	ALIGN(r5)
	ALIGN(r6)
	ALIGN(r7)
	ALIGN(r8)
	ALIGN(r9)
	ALIGN(r10)
	ALIGN(r11)
	ALIGN(r12)
	ALIGN(r13)
	ALIGN(r14)
	ALIGN(r15)
	
	UNPACK256_LO_HI_EPI16(0, 1)
	UNPACK256_LO_HI_EPI16(2, 3)
	UNPACK256_LO_HI_EPI16(4, 5)
	UNPACK256_LO_HI_EPI16(6, 7)
	UNPACK256_LO_HI_EPI16(8, 9)
	UNPACK256_LO_HI_EPI16(10, 11)
	UNPACK256_LO_HI_EPI16(12, 13)
	UNPACK256_LO_HI_EPI16(14, 15)

	UNPACK256_LO_HI_EPI32(0, 2)
	UNPACK256_LO_HI_EPI32(1, 3)
	UNPACK256_LO_HI_EPI32(4, 6)
	UNPACK256_LO_HI_EPI32(5, 7)
	UNPACK256_LO_HI_EPI32(8, 10)
	UNPACK256_LO_HI_EPI32(9, 11)
	UNPACK256_LO_HI_EPI32(12, 14)
	UNPACK256_LO_HI_EPI32(13, 15)

	UNPACK256_LO_HI_EPI64(0, 4)
	UNPACK256_LO_HI_EPI64(2, 6)
	UNPACK256_LO_HI_EPI64(1, 5)
	UNPACK256_LO_HI_EPI64(3, 7)
	UNPACK256_LO_HI_EPI64(8, 12)
	UNPACK256_LO_HI_EPI64(10, 14)
	UNPACK256_LO_HI_EPI64(9, 13)
	UNPACK256_LO_HI_EPI64(11, 15)

	UNPACK256_LO_HI_EPI128(0, 8)
	UNPACK256_LO_HI_EPI128(4, 12)
	UNPACK256_LO_HI_EPI128(2, 10)
	UNPACK256_LO_HI_EPI128(6, 14)
	UNPACK256_LO_HI_EPI128(1, 9)
	UNPACK256_LO_HI_EPI128(5, 13)
	UNPACK256_LO_HI_EPI128(3, 11)
	UNPACK256_LO_HI_EPI128(7, 15)
	
	__m256i* ptr = (__m256i*)out;
	_mm256_store_si256(ptr++, r0);
	_mm256_store_si256(ptr++, r4);
	_mm256_store_si256(ptr++, r2);
	_mm256_store_si256(ptr++, r6);
	_mm256_store_si256(ptr++, r1);
	_mm256_store_si256(ptr++, r5);
	_mm256_store_si256(ptr++, r3);
	_mm256_store_si256(ptr++, r7);
	_mm256_store_si256(ptr++, r8);
	_mm256_store_si256(ptr++, r12);
	_mm256_store_si256(ptr++, r10);
	_mm256_store_si256(ptr++, r14);
	_mm256_store_si256(ptr++, r9);
	_mm256_store_si256(ptr++, r13);
	_mm256_store_si256(ptr++, r11);
	_mm256_store_si256(ptr++, r15);
}

#endif

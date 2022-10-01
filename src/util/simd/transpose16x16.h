/****
DIAMOND protein aligner
Copyright (C) 2020 Max Planck Society for the Advancement of Science e.V.

Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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

#pragma once
#include "../simd.h"

#if defined(__SSE2__)

#define UNPACK128_LO_HI_EPI8(a, b) t = r##a; r##a = _mm_unpacklo_epi8(t, r##b); r##b = _mm_unpackhi_epi8(t, r##b);
#define UNPACK128_LO_HI_EPI16(a, b) t = r##a; r##a = _mm_unpacklo_epi16(t, r##b); r##b = _mm_unpackhi_epi16(t, r##b);
#define UNPACK128_LO_HI_EPI32(a, b) t = r##a; r##a = _mm_unpacklo_epi32(t, r##b); r##b = _mm_unpackhi_epi32(t, r##b);
#define UNPACK128_LO_HI_EPI64(a, b) t = r##a; r##a = _mm_unpacklo_epi64(t, r##b); r##b = _mm_unpackhi_epi64(t, r##b);

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
#define UNPACK128_LO_HI_EPI64(a, b) t0 = vget_high_s8(r##a); \
        t1 = vget_low_s8(r##a); \
        t2 = vget_high_s8(r##b); \
        t3 = vget_low_s8(r##b); \
        r##a = vcombine_s8(t1, t3); \
        r##b = vcombine_s8(t0, t2);

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


	b0.val[0] = vcombine_s8(vget_low_s8(r0), vget_low_s8(r8));
	b0.val[1] = vcombine_s8(vget_low_s8(r1), vget_low_s8(r9));
	b0.val[2] = vcombine_s8(vget_low_s8(r2), vget_low_s8(r10));
	b0.val[3] = vcombine_s8(vget_low_s8(r3), vget_low_s8(r11));
	vst1q_s8_x4(&out[0x00], b0);

	b1.val[0] = vcombine_s8(vget_low_s8(r4), vget_low_s8(r12));
	b1.val[1] = vcombine_s8(vget_low_s8(r5), vget_low_s8(r13));
	b1.val[2] = vcombine_s8(vget_low_s8(r6), vget_low_s8(r14));
	b1.val[3] = vcombine_s8(vget_low_s8(r7), vget_low_s8(r15));
	vst1q_s8_x4(&out[0x40], b1);

	b2.val[0] = vcombine_s8(vget_high_s8(r0), vget_high_s8(r8));
	b2.val[1] = vcombine_s8(vget_high_s8(r1), vget_high_s8(r9));
	b2.val[2] = vcombine_s8(vget_high_s8(r2), vget_high_s8(r10));
	b2.val[3] = vcombine_s8(vget_high_s8(r3), vget_high_s8(r11));
        vst1q_s8_x4(&out[0x80], b2);

	b3.val[0] = vcombine_s8(vget_high_s8(r4), vget_high_s8(r12));
	b3.val[1] = vcombine_s8(vget_high_s8(r5), vget_high_s8(r13));
	b3.val[2] = vcombine_s8(vget_high_s8(r6), vget_high_s8(r14));
	b3.val[3] = vcombine_s8(vget_high_s8(r7), vget_high_s8(r15));
	vst1q_s8_x4(&out[0xc0], b3);
}

#endif

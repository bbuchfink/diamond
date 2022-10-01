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

#elif defined(__aarch64__)

#define UNPACK128_LO_HI_EPI8(a, b) t = r##a; \
	r##a = vzip1q_s8(t, r##b); \
	r##b = vzip2q_s8(t, r##b);
#define UNPACK128_LO_HI_EPI16(a, b) t = r##a; \
	r##a = vreinterpretq_s8_s16(vzip1q_s16(vreinterpretq_s16_s8(t), vreinterpretq_s16_s8(r##b))); \
	r##b = vreinterpretq_s8_s16(vzip2q_s16(vreinterpretq_s16_s8(t), vreinterpretq_s16_s8(r##b)));
#define UNPACK128_LO_HI_EPI32(a, b) t = r##a; \
	r##a = vreinterpretq_s8_s32(vzip1q_s32(vreinterpretq_s32_s8(t), vreinterpretq_s32_s8(r##b))); \
	r##b = vreinterpretq_s8_s32(vzip2q_s32(vreinterpretq_s32_s8(t), vreinterpretq_s32_s8(r##b)));
#define UNPACK128_LO_HI_EPI64(a, b) t = r##a; \
	r##a = vreinterpretq_s8_s64(vzip1q_s64(vreinterpretq_s64_s8(t), vreinterpretq_s64_s8(r##b))); \
	r##b = vreinterpretq_s8_s64(vzip2q_s64(vreinterpretq_s64_s8(t), vreinterpretq_s64_s8(r##b)));

static inline void transpose(const signed char **data, size_t n, signed char *out, const int8x16_t&) {
	int8x16_t r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15, t;
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

	UNPACK128_LO_HI_EPI64(0, 8)
	UNPACK128_LO_HI_EPI64(4, 12)
	UNPACK128_LO_HI_EPI64(2, 10)
	UNPACK128_LO_HI_EPI64(6, 14)
	UNPACK128_LO_HI_EPI64(1, 9)
	UNPACK128_LO_HI_EPI64(5, 13)
	UNPACK128_LO_HI_EPI64(3, 11)
	UNPACK128_LO_HI_EPI64(7, 15)

	int8x16_t* ptr = (int8x16_t*)out;
	vst1q_s8((int8_t*) (ptr++), r0 );
	vst1q_s8((int8_t*) (ptr++), r8 );
	vst1q_s8((int8_t*) (ptr++), r4 );
	vst1q_s8((int8_t*) (ptr++), r12);
	vst1q_s8((int8_t*) (ptr++), r2 );
	vst1q_s8((int8_t*) (ptr++), r10);
	vst1q_s8((int8_t*) (ptr++), r6 );
	vst1q_s8((int8_t*) (ptr++), r14);
	vst1q_s8((int8_t*) (ptr++), r1 );
	vst1q_s8((int8_t*) (ptr++), r9 );
	vst1q_s8((int8_t*) (ptr++), r5 );
	vst1q_s8((int8_t*) (ptr++), r13);
	vst1q_s8((int8_t*) (ptr++), r3 );
	vst1q_s8((int8_t*) (ptr++), r11);
	vst1q_s8((int8_t*) (ptr++), r7 );
	vst1q_s8((int8_t*) (ptr),   r15);
}

#endif

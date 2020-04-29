/****
DIAMOND protein aligner
Copyright (C) 2020 Max Planck Society for the Advancement of Science e.V.

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

#ifndef TRANSPOSE16X16_H_
#define TRANSPOSE16X16_H_

#include "../simd.h"

static inline void transpose16x16(const signed char **data, size_t n, signed char *out) {
	__m128i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15;
	
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

	__m128i t = r0;
	r0 = _mm_unpacklo_epi8(t, r1);
	r1 = _mm_unpackhi_epi8(t, r1);
	t = r2;
	r2 = _mm_unpacklo_epi8(t, r3);
	r3 = _mm_unpackhi_epi8(t, r3);
	t = r4;
	r4 = _mm_unpacklo_epi8(t, r5);
	r5 = _mm_unpackhi_epi8(t, r5);
	t = r6;
	r6 = _mm_unpacklo_epi8(t, r7);
	r7 = _mm_unpackhi_epi8(t, r7);
	t = r8;
	r8 = _mm_unpacklo_epi8(t, r9);
	r9 = _mm_unpackhi_epi8(t, r9);
	t = r10;
	r10 = _mm_unpacklo_epi8(t, r11);
	r11 = _mm_unpackhi_epi8(t, r11);
	t = r12;
	r12 = _mm_unpacklo_epi8(t, r13);
	r13 = _mm_unpackhi_epi8(t, r13);
	t = r14;
	r14 = _mm_unpacklo_epi8(t, r15);
	r15 = _mm_unpackhi_epi8(t, r15);

	t = r0;
	r0 = _mm_unpacklo_epi16(t, r2);
	r2 = _mm_unpackhi_epi16(t, r2);
	t = r4;
	r4 = _mm_unpacklo_epi16(t, r6);
	r6 = _mm_unpackhi_epi16(t, r6);
	t = r8;
	r8 = _mm_unpacklo_epi16(t, r10);
	r10 = _mm_unpackhi_epi16(t, r10);
	t = r12;
	r12 = _mm_unpacklo_epi16(t, r14);
	r14 = _mm_unpackhi_epi16(t, r14);
	t = r1;
	r1 = _mm_unpacklo_epi16(t, r3);
	r3 = _mm_unpackhi_epi16(t, r3);
	t = r5;
	r5 = _mm_unpacklo_epi16(t, r7);
	r7 = _mm_unpackhi_epi16(t, r7);
	t = r9;
	r9 = _mm_unpacklo_epi16(t, r11);
	r11 = _mm_unpackhi_epi16(t, r11);
	t = r13;
	r13 = _mm_unpacklo_epi16(t, r15);
	r15 = _mm_unpackhi_epi16(t, r15);

	t = r0;
	r0 = _mm_unpacklo_epi32(t, r4);
	r4 = _mm_unpackhi_epi32(t, r4);
	t = r8;
	r8 = _mm_unpacklo_epi32(t, r12);
	r12 = _mm_unpackhi_epi32(t, r12);
	t = r2;
	r2 = _mm_unpacklo_epi32(t, r6);
	r6 = _mm_unpackhi_epi32(t, r6);
	t = r10;
	r10 = _mm_unpacklo_epi32(t, r14);
	r14 = _mm_unpackhi_epi32(t, r14);
	t = r1;
	r1 = _mm_unpacklo_epi32(t, r5);
	r5 = _mm_unpackhi_epi32(t, r5);
	t = r9;
	r9 = _mm_unpacklo_epi32(t, r13);
	r13 = _mm_unpackhi_epi32(t, r13);
	t = r3;
	r3 = _mm_unpacklo_epi32(t, r7);
	r7 = _mm_unpackhi_epi32(t, r7);
	t = r11;
	r11 = _mm_unpacklo_epi32(t, r15);
	r15 = _mm_unpackhi_epi32(t, r15);

	t = r0;
	r0 = _mm_unpacklo_epi64(t, r8);
	r8 = _mm_unpackhi_epi64(t, r8);
	t = r4;
	r4 = _mm_unpacklo_epi64(t, r12);
	r12 = _mm_unpackhi_epi64(t, r12);
	t = r2;
	r2 = _mm_unpacklo_epi64(t, r10);
	r10 = _mm_unpackhi_epi64(t, r10);
	t = r6;
	r6 = _mm_unpacklo_epi64(t, r14);
	r14 = _mm_unpackhi_epi64(t, r14);
	t = r1;
	r1 = _mm_unpacklo_epi64(t, r9);
	r9 = _mm_unpackhi_epi64(t, r9);
	t = r5;
	r5 = _mm_unpacklo_epi64(t, r13);
	r13 = _mm_unpackhi_epi64(t, r13);
	t = r3;
	r3 = _mm_unpacklo_epi64(t, r11);
	r11 = _mm_unpackhi_epi64(t, r11);
	t = r7;
	r7 = _mm_unpacklo_epi64(t, r15);
	r15 = _mm_unpackhi_epi64(t, r15);

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

#endif
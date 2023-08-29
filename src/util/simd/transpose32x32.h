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
#include <stddef.h>
#include "../simd.h"

#define UNPACK256_LO_HI_EPI8(a, b) t = r##a; r##a = _mm256_unpacklo_epi8(t, r##b); r##b = _mm256_unpackhi_epi8(t, r##b);
#define UNPACK256_LO_HI_EPI16(a, b) t = r##a; r##a = _mm256_unpacklo_epi16(t, r##b); r##b = _mm256_unpackhi_epi16(t, r##b);
#define UNPACK256_LO_HI_EPI32(a, b) t = r##a; r##a = _mm256_unpacklo_epi32(t, r##b); r##b = _mm256_unpackhi_epi32(t, r##b);
#define UNPACK256_LO_HI_EPI64(a, b) t = r##a; r##a = _mm256_unpacklo_epi64(t, r##b); r##b = _mm256_unpackhi_epi64(t, r##b);
#define UNPACK256_LO_HI_EPI128(a, b) t = r##a; r##a = _mm256_permute2x128_si256(t, r##b, 32); r##b = _mm256_permute2x128_si256(t, r##b, 49);

static inline void transpose(const signed char** data, size_t n, signed char* out, const __m256i&) {
	__m256i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15, r16, r17, r18, r19, r20, r21, r22, r23, r24, r25, r26, r27, r28, r29, r30, r31, t;
	r0 = r1 = r2 = r3 = r4 = r5 = r6 = r7 = r8 = r9 = r10 = r11 = r12 = r13 = r14 = r15 = r16 = r17 = r18 = r19 = r20 = r21 = r22 = r23 = r24 = r25 = r26 = r27 = r28 = r29 = r30 = r31 = _mm256_setzero_si256();

	switch (n) {
	case 32: r0 = _mm256_loadu_si256((const __m256i*) * (data++));
	case 31: r1 = _mm256_loadu_si256((const __m256i*) * (data++));
	case 30: r2 = _mm256_loadu_si256((const __m256i*) * (data++));
	case 29: r3 = _mm256_loadu_si256((const __m256i*) * (data++));
	case 28: r4 = _mm256_loadu_si256((const __m256i*) * (data++));
	case 27: r5 = _mm256_loadu_si256((const __m256i*) * (data++));
	case 26: r6 = _mm256_loadu_si256((const __m256i*) * (data++));
	case 25: r7 = _mm256_loadu_si256((const __m256i*) * (data++));
	case 24: r8 = _mm256_loadu_si256((const __m256i*) * (data++));
	case 23: r9 = _mm256_loadu_si256((const __m256i*) * (data++));
	case 22: r10 = _mm256_loadu_si256((const __m256i*) * (data++));
	case 21: r11 = _mm256_loadu_si256((const __m256i*) * (data++));
	case 20: r12 = _mm256_loadu_si256((const __m256i*) * (data++));
	case 19: r13 = _mm256_loadu_si256((const __m256i*) * (data++));
	case 18: r14 = _mm256_loadu_si256((const __m256i*) * (data++));
	case 17: r15 = _mm256_loadu_si256((const __m256i*) * (data++));
	case 16: r16 = _mm256_loadu_si256((const __m256i*) * (data++));
	case 15: r17 = _mm256_loadu_si256((const __m256i*) * (data++));
	case 14: r18 = _mm256_loadu_si256((const __m256i*) * (data++));
	case 13: r19 = _mm256_loadu_si256((const __m256i*) * (data++));
	case 12: r20 = _mm256_loadu_si256((const __m256i*) * (data++));
	case 11: r21 = _mm256_loadu_si256((const __m256i*) * (data++));
	case 10: r22 = _mm256_loadu_si256((const __m256i*) * (data++));
	case 9: r23 = _mm256_loadu_si256((const __m256i*) * (data++));
	case 8: r24 = _mm256_loadu_si256((const __m256i*) * (data++));
	case 7: r25 = _mm256_loadu_si256((const __m256i*) * (data++));
	case 6: r26 = _mm256_loadu_si256((const __m256i*) * (data++));
	case 5: r27 = _mm256_loadu_si256((const __m256i*) * (data++));
	case 4: r28 = _mm256_loadu_si256((const __m256i*) * (data++));
	case 3: r29 = _mm256_loadu_si256((const __m256i*) * (data++));
	case 2: r30 = _mm256_loadu_si256((const __m256i*) * (data++));
	case 1: r31 = _mm256_loadu_si256((const __m256i*) * data);
	}
	
	UNPACK256_LO_HI_EPI8(0, 1)
	UNPACK256_LO_HI_EPI8(2, 3)
	UNPACK256_LO_HI_EPI8(4, 5)
	UNPACK256_LO_HI_EPI8(6, 7)
	UNPACK256_LO_HI_EPI8(8, 9)
	UNPACK256_LO_HI_EPI8(10, 11)
	UNPACK256_LO_HI_EPI8(12, 13)
	UNPACK256_LO_HI_EPI8(14, 15)
	UNPACK256_LO_HI_EPI8(16, 17)
	UNPACK256_LO_HI_EPI8(18, 19)
	UNPACK256_LO_HI_EPI8(20, 21)
	UNPACK256_LO_HI_EPI8(22, 23)
	UNPACK256_LO_HI_EPI8(24, 25)
	UNPACK256_LO_HI_EPI8(26, 27)
	UNPACK256_LO_HI_EPI8(28, 29)
	UNPACK256_LO_HI_EPI8(30, 31)

	UNPACK256_LO_HI_EPI16(0, 2)
	UNPACK256_LO_HI_EPI16(1, 3)
	UNPACK256_LO_HI_EPI16(4, 6)
	UNPACK256_LO_HI_EPI16(5, 7)
	UNPACK256_LO_HI_EPI16(8, 10)
	UNPACK256_LO_HI_EPI16(9, 11)
	UNPACK256_LO_HI_EPI16(12, 14)
	UNPACK256_LO_HI_EPI16(13, 15)
	UNPACK256_LO_HI_EPI16(16, 18)
	UNPACK256_LO_HI_EPI16(17, 19)
	UNPACK256_LO_HI_EPI16(20, 22)
	UNPACK256_LO_HI_EPI16(21, 23)
	UNPACK256_LO_HI_EPI16(24, 26)
	UNPACK256_LO_HI_EPI16(25, 27)
	UNPACK256_LO_HI_EPI16(28, 30)
	UNPACK256_LO_HI_EPI16(29, 31)

	UNPACK256_LO_HI_EPI32(0, 4)
	UNPACK256_LO_HI_EPI32(2, 6)
	UNPACK256_LO_HI_EPI32(1, 5)
	UNPACK256_LO_HI_EPI32(3, 7)
	UNPACK256_LO_HI_EPI32(8, 12)
	UNPACK256_LO_HI_EPI32(10, 14)
	UNPACK256_LO_HI_EPI32(9, 13)
	UNPACK256_LO_HI_EPI32(11, 15)
	UNPACK256_LO_HI_EPI32(16, 20)
	UNPACK256_LO_HI_EPI32(18, 22)
	UNPACK256_LO_HI_EPI32(17, 21)
	UNPACK256_LO_HI_EPI32(19, 23)
	UNPACK256_LO_HI_EPI32(24, 28)
	UNPACK256_LO_HI_EPI32(26, 30)
	UNPACK256_LO_HI_EPI32(25, 29)
	UNPACK256_LO_HI_EPI32(27, 31)

	UNPACK256_LO_HI_EPI64(0, 8)
	UNPACK256_LO_HI_EPI64(4, 12)
	UNPACK256_LO_HI_EPI64(2, 10)
	UNPACK256_LO_HI_EPI64(6, 14)
	UNPACK256_LO_HI_EPI64(1, 9)
	UNPACK256_LO_HI_EPI64(5, 13)
	UNPACK256_LO_HI_EPI64(3, 11)
	UNPACK256_LO_HI_EPI64(7, 15)
	UNPACK256_LO_HI_EPI64(16, 24)
	UNPACK256_LO_HI_EPI64(20, 28)
	UNPACK256_LO_HI_EPI64(18, 26)
	UNPACK256_LO_HI_EPI64(22, 30)
	UNPACK256_LO_HI_EPI64(17, 25)
	UNPACK256_LO_HI_EPI64(21, 29)
	UNPACK256_LO_HI_EPI64(19, 27)
	UNPACK256_LO_HI_EPI64(23, 31)

	UNPACK256_LO_HI_EPI128(0, 16)
	UNPACK256_LO_HI_EPI128(8, 24)
	UNPACK256_LO_HI_EPI128(4, 20)
	UNPACK256_LO_HI_EPI128(12, 28)
	UNPACK256_LO_HI_EPI128(2, 18)
	UNPACK256_LO_HI_EPI128(10, 26)
	UNPACK256_LO_HI_EPI128(6, 22)
	UNPACK256_LO_HI_EPI128(14, 30)
	UNPACK256_LO_HI_EPI128(1, 17)
	UNPACK256_LO_HI_EPI128(9, 25)
	UNPACK256_LO_HI_EPI128(5, 21)
	UNPACK256_LO_HI_EPI128(13, 29)
	UNPACK256_LO_HI_EPI128(3, 19)
	UNPACK256_LO_HI_EPI128(11, 27)
	UNPACK256_LO_HI_EPI128(7, 23)
	UNPACK256_LO_HI_EPI128(15, 31)
	
	__m256i* ptr = (__m256i*)out;
	_mm256_store_si256(ptr++, r0);
	_mm256_store_si256(ptr++, r8);
	_mm256_store_si256(ptr++, r4);
	_mm256_store_si256(ptr++, r12);
	_mm256_store_si256(ptr++, r2);
	_mm256_store_si256(ptr++, r10);
	_mm256_store_si256(ptr++, r6);
	_mm256_store_si256(ptr++, r14);
	_mm256_store_si256(ptr++, r1);
	_mm256_store_si256(ptr++, r9);
	_mm256_store_si256(ptr++, r5);
	_mm256_store_si256(ptr++, r13);
	_mm256_store_si256(ptr++, r3);
	_mm256_store_si256(ptr++, r11);
	_mm256_store_si256(ptr++, r7);
	_mm256_store_si256(ptr++, r15);
	_mm256_store_si256(ptr++, r16);
	_mm256_store_si256(ptr++, r24);
	_mm256_store_si256(ptr++, r20);
	_mm256_store_si256(ptr++, r28);
	_mm256_store_si256(ptr++, r18);
	_mm256_store_si256(ptr++, r26);
	_mm256_store_si256(ptr++, r22);
	_mm256_store_si256(ptr++, r30);
	_mm256_store_si256(ptr++, r17);
	_mm256_store_si256(ptr++, r25);
	_mm256_store_si256(ptr++, r21);
	_mm256_store_si256(ptr++, r29);
	_mm256_store_si256(ptr++, r19);
	_mm256_store_si256(ptr++, r27);
	_mm256_store_si256(ptr++, r23);
	_mm256_store_si256(ptr++, r31);
}

static inline void transpose_offset(const signed char** data, size_t n, ptrdiff_t offset, signed char* out, __m256i) {
	__m256i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15, r16, r17, r18, r19, r20, r21, r22, r23, r24, r25, r26, r27, r28, r29, r30, r31, t;
	r0 = r1 = r2 = r3 = r4 = r5 = r6 = r7 = r8 = r9 = r10 = r11 = r12 = r13 = r14 = r15 = r16 = r17 = r18 = r19 = r20 = r21 = r22 = r23 = r24 = r25 = r26 = r27 = r28 = r29 = r30 = r31 = _mm256_setzero_si256();

	switch (n) {
	case 32: r0 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	case 31: r1 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	case 30: r2 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	case 29: r3 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	case 28: r4 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	case 27: r5 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	case 26: r6 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	case 25: r7 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	case 24: r8 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	case 23: r9 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	case 22: r10 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	case 21: r11 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	case 20: r12 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	case 19: r13 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	case 18: r14 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	case 17: r15 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	case 16: r16 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	case 15: r17 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	case 14: r18 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	case 13: r19 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	case 12: r20 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	case 11: r21 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	case 10: r22 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	case 9: r23 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	case 8: r24 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	case 7: r25 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	case 6: r26 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	case 5: r27 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	case 4: r28 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	case 3: r29 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	case 2: r30 = _mm256_loadu_si256((const __m256i*) * (data++) + offset);
	case 1: r31 = _mm256_loadu_si256((const __m256i*) * data + offset);
	}
	
	UNPACK256_LO_HI_EPI8(0, 1)
	UNPACK256_LO_HI_EPI8(2, 3)
	UNPACK256_LO_HI_EPI8(4, 5)
	UNPACK256_LO_HI_EPI8(6, 7)
	UNPACK256_LO_HI_EPI8(8, 9)
	UNPACK256_LO_HI_EPI8(10, 11)
	UNPACK256_LO_HI_EPI8(12, 13)
	UNPACK256_LO_HI_EPI8(14, 15)
	UNPACK256_LO_HI_EPI8(16, 17)
	UNPACK256_LO_HI_EPI8(18, 19)
	UNPACK256_LO_HI_EPI8(20, 21)
	UNPACK256_LO_HI_EPI8(22, 23)
	UNPACK256_LO_HI_EPI8(24, 25)
	UNPACK256_LO_HI_EPI8(26, 27)
	UNPACK256_LO_HI_EPI8(28, 29)
	UNPACK256_LO_HI_EPI8(30, 31)

	UNPACK256_LO_HI_EPI16(0, 2)
	UNPACK256_LO_HI_EPI16(1, 3)
	UNPACK256_LO_HI_EPI16(4, 6)
	UNPACK256_LO_HI_EPI16(5, 7)
	UNPACK256_LO_HI_EPI16(8, 10)
	UNPACK256_LO_HI_EPI16(9, 11)
	UNPACK256_LO_HI_EPI16(12, 14)
	UNPACK256_LO_HI_EPI16(13, 15)
	UNPACK256_LO_HI_EPI16(16, 18)
	UNPACK256_LO_HI_EPI16(17, 19)
	UNPACK256_LO_HI_EPI16(20, 22)
	UNPACK256_LO_HI_EPI16(21, 23)
	UNPACK256_LO_HI_EPI16(24, 26)
	UNPACK256_LO_HI_EPI16(25, 27)
	UNPACK256_LO_HI_EPI16(28, 30)
	UNPACK256_LO_HI_EPI16(29, 31)

	UNPACK256_LO_HI_EPI32(0, 4)
	UNPACK256_LO_HI_EPI32(2, 6)
	UNPACK256_LO_HI_EPI32(1, 5)
	UNPACK256_LO_HI_EPI32(3, 7)
	UNPACK256_LO_HI_EPI32(8, 12)
	UNPACK256_LO_HI_EPI32(10, 14)
	UNPACK256_LO_HI_EPI32(9, 13)
	UNPACK256_LO_HI_EPI32(11, 15)
	UNPACK256_LO_HI_EPI32(16, 20)
	UNPACK256_LO_HI_EPI32(18, 22)
	UNPACK256_LO_HI_EPI32(17, 21)
	UNPACK256_LO_HI_EPI32(19, 23)
	UNPACK256_LO_HI_EPI32(24, 28)
	UNPACK256_LO_HI_EPI32(26, 30)
	UNPACK256_LO_HI_EPI32(25, 29)
	UNPACK256_LO_HI_EPI32(27, 31)

	UNPACK256_LO_HI_EPI64(0, 8)
	UNPACK256_LO_HI_EPI64(4, 12)
	UNPACK256_LO_HI_EPI64(2, 10)
	UNPACK256_LO_HI_EPI64(6, 14)
	UNPACK256_LO_HI_EPI64(1, 9)
	UNPACK256_LO_HI_EPI64(5, 13)
	UNPACK256_LO_HI_EPI64(3, 11)
	UNPACK256_LO_HI_EPI64(7, 15)
	UNPACK256_LO_HI_EPI64(16, 24)
	UNPACK256_LO_HI_EPI64(20, 28)
	UNPACK256_LO_HI_EPI64(18, 26)
	UNPACK256_LO_HI_EPI64(22, 30)
	UNPACK256_LO_HI_EPI64(17, 25)
	UNPACK256_LO_HI_EPI64(21, 29)
	UNPACK256_LO_HI_EPI64(19, 27)
	UNPACK256_LO_HI_EPI64(23, 31)

	UNPACK256_LO_HI_EPI128(0, 16)
	UNPACK256_LO_HI_EPI128(8, 24)
	UNPACK256_LO_HI_EPI128(4, 20)
	UNPACK256_LO_HI_EPI128(12, 28)
	UNPACK256_LO_HI_EPI128(2, 18)
	UNPACK256_LO_HI_EPI128(10, 26)
	UNPACK256_LO_HI_EPI128(6, 22)
	UNPACK256_LO_HI_EPI128(14, 30)
	UNPACK256_LO_HI_EPI128(1, 17)
	UNPACK256_LO_HI_EPI128(9, 25)
	UNPACK256_LO_HI_EPI128(5, 21)
	UNPACK256_LO_HI_EPI128(13, 29)
	UNPACK256_LO_HI_EPI128(3, 19)
	UNPACK256_LO_HI_EPI128(11, 27)
	UNPACK256_LO_HI_EPI128(7, 23)
	UNPACK256_LO_HI_EPI128(15, 31)
	
	__m256i* ptr = (__m256i*)out;
	_mm256_store_si256(ptr++, r0);
	_mm256_store_si256(ptr++, r8);
	_mm256_store_si256(ptr++, r4);
	_mm256_store_si256(ptr++, r12);
	_mm256_store_si256(ptr++, r2);
	_mm256_store_si256(ptr++, r10);
	_mm256_store_si256(ptr++, r6);
	_mm256_store_si256(ptr++, r14);
	_mm256_store_si256(ptr++, r1);
	_mm256_store_si256(ptr++, r9);
	_mm256_store_si256(ptr++, r5);
	_mm256_store_si256(ptr++, r13);
	_mm256_store_si256(ptr++, r3);
	_mm256_store_si256(ptr++, r11);
	_mm256_store_si256(ptr++, r7);
	_mm256_store_si256(ptr++, r15);
	_mm256_store_si256(ptr++, r16);
	_mm256_store_si256(ptr++, r24);
	_mm256_store_si256(ptr++, r20);
	_mm256_store_si256(ptr++, r28);
	_mm256_store_si256(ptr++, r18);
	_mm256_store_si256(ptr++, r26);
	_mm256_store_si256(ptr++, r22);
	_mm256_store_si256(ptr++, r30);
	_mm256_store_si256(ptr++, r17);
	_mm256_store_si256(ptr++, r25);
	_mm256_store_si256(ptr++, r21);
	_mm256_store_si256(ptr++, r29);
	_mm256_store_si256(ptr++, r19);
	_mm256_store_si256(ptr++, r27);
	_mm256_store_si256(ptr++, r23);
	_mm256_store_si256(ptr++, r31);
}
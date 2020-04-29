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

#ifndef TRANSPOSE32X32_H_
#define TRANSPOSE32X32_H_

#include "../simd.h"

static inline void transpose32x32(const signed char** data, size_t n, signed char* out) {
	__m256i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15, r16, r17, r18, r19, r20, r21, r22, r23, r24, r25, r26, r27, r28, r29, r30, r31;

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

	__m256i t = r0;
	r0 = _mm256_unpacklo_epi8(t, r1);
	r1 = _mm256_unpackhi_epi8(t, r1);
	t = r2;
	r2 = _mm256_unpacklo_epi8(t, r3);
	r3 = _mm256_unpackhi_epi8(t, r3);
	t = r4;
	r4 = _mm256_unpacklo_epi8(t, r5);
	r5 = _mm256_unpackhi_epi8(t, r5);
	t = r6;
	r6 = _mm256_unpacklo_epi8(t, r7);
	r7 = _mm256_unpackhi_epi8(t, r7);
	t = r8;
	r8 = _mm256_unpacklo_epi8(t, r9);
	r9 = _mm256_unpackhi_epi8(t, r9);
	t = r10;
	r10 = _mm256_unpacklo_epi8(t, r11);
	r11 = _mm256_unpackhi_epi8(t, r11);
	t = r12;
	r12 = _mm256_unpacklo_epi8(t, r13);
	r13 = _mm256_unpackhi_epi8(t, r13);
	t = r14;
	r14 = _mm256_unpacklo_epi8(t, r15);
	r15 = _mm256_unpackhi_epi8(t, r15);
	t = r16;
	r16 = _mm256_unpacklo_epi8(t, r17);
	r17 = _mm256_unpackhi_epi8(t, r17);
	t = r18;
	r18 = _mm256_unpacklo_epi8(t, r19);
	r19 = _mm256_unpackhi_epi8(t, r19);
	t = r20;
	r20 = _mm256_unpacklo_epi8(t, r21);
	r21 = _mm256_unpackhi_epi8(t, r21);
	t = r22;
	r22 = _mm256_unpacklo_epi8(t, r23);
	r23 = _mm256_unpackhi_epi8(t, r23);
	t = r24;
	r24 = _mm256_unpacklo_epi8(t, r25);
	r25 = _mm256_unpackhi_epi8(t, r25);
	t = r26;
	r26 = _mm256_unpacklo_epi8(t, r27);
	r27 = _mm256_unpackhi_epi8(t, r27);
	t = r28;
	r28 = _mm256_unpacklo_epi8(t, r29);
	r29 = _mm256_unpackhi_epi8(t, r29);
	t = r30;
	r30 = _mm256_unpacklo_epi8(t, r31);
	r31 = _mm256_unpackhi_epi8(t, r31);

	t = r0;
	r0 = _mm256_unpacklo_epi16(t, r2);
	r2 = _mm256_unpackhi_epi16(t, r2);
	t = r4;
	r4 = _mm256_unpacklo_epi16(t, r6);
	r6 = _mm256_unpackhi_epi16(t, r6);
	t = r8;
	r8 = _mm256_unpacklo_epi16(t, r10);
	r10 = _mm256_unpackhi_epi16(t, r10);
	t = r12;
	r12 = _mm256_unpacklo_epi16(t, r14);
	r14 = _mm256_unpackhi_epi16(t, r14);
	t = r16;
	r16 = _mm256_unpacklo_epi16(t, r18);
	r18 = _mm256_unpackhi_epi16(t, r18);
	t = r20;
	r20 = _mm256_unpacklo_epi16(t, r22);
	r22 = _mm256_unpackhi_epi16(t, r22);
	t = r24;
	r24 = _mm256_unpacklo_epi16(t, r26);
	r26 = _mm256_unpackhi_epi16(t, r26);
	t = r28;
	r28 = _mm256_unpacklo_epi16(t, r30);
	r30 = _mm256_unpackhi_epi16(t, r30);
	t = r1;
	r1 = _mm256_unpacklo_epi16(t, r3);
	r3 = _mm256_unpackhi_epi16(t, r3);
	t = r5;
	r5 = _mm256_unpacklo_epi16(t, r7);
	r7 = _mm256_unpackhi_epi16(t, r7);
	t = r9;
	r9 = _mm256_unpacklo_epi16(t, r11);
	r11 = _mm256_unpackhi_epi16(t, r11);
	t = r13;
	r13 = _mm256_unpacklo_epi16(t, r15);
	r15 = _mm256_unpackhi_epi16(t, r15);
	t = r17;
	r17 = _mm256_unpacklo_epi16(t, r19);
	r19 = _mm256_unpackhi_epi16(t, r19);
	t = r21;
	r21 = _mm256_unpacklo_epi16(t, r23);
	r23 = _mm256_unpackhi_epi16(t, r23);
	t = r25;
	r25 = _mm256_unpacklo_epi16(t, r27);
	r27 = _mm256_unpackhi_epi16(t, r27);
	t = r29;
	r29 = _mm256_unpacklo_epi16(t, r31);
	r31 = _mm256_unpackhi_epi16(t, r31);

	t = r0;
	r0 = _mm256_unpacklo_epi32(t, r4);
	r4 = _mm256_unpackhi_epi32(t, r4);
	t = r8;
	r8 = _mm256_unpacklo_epi32(t, r12);
	r12 = _mm256_unpackhi_epi32(t, r12);
	t = r16;
	r16 = _mm256_unpacklo_epi32(t, r20);
	r20 = _mm256_unpackhi_epi32(t, r20);
	t = r24;
	r24 = _mm256_unpacklo_epi32(t, r28);
	r28 = _mm256_unpackhi_epi32(t, r28);
	t = r2;
	r2 = _mm256_unpacklo_epi32(t, r6);
	r6 = _mm256_unpackhi_epi32(t, r6);
	t = r10;
	r10 = _mm256_unpacklo_epi32(t, r14);
	r14 = _mm256_unpackhi_epi32(t, r14);
	t = r18;
	r18 = _mm256_unpacklo_epi32(t, r22);
	r22 = _mm256_unpackhi_epi32(t, r22);
	t = r26;
	r26 = _mm256_unpacklo_epi32(t, r30);
	r30 = _mm256_unpackhi_epi32(t, r30);
	t = r1;
	r1 = _mm256_unpacklo_epi32(t, r5);
	r5 = _mm256_unpackhi_epi32(t, r5);
	t = r9;
	r9 = _mm256_unpacklo_epi32(t, r13);
	r13 = _mm256_unpackhi_epi32(t, r13);
	t = r17;
	r17 = _mm256_unpacklo_epi32(t, r21);
	r21 = _mm256_unpackhi_epi32(t, r21);
	t = r25;
	r25 = _mm256_unpacklo_epi32(t, r29);
	r29 = _mm256_unpackhi_epi32(t, r29);
	t = r3;
	r3 = _mm256_unpacklo_epi32(t, r7);
	r7 = _mm256_unpackhi_epi32(t, r7);
	t = r11;
	r11 = _mm256_unpacklo_epi32(t, r15);
	r15 = _mm256_unpackhi_epi32(t, r15);
	t = r19;
	r19 = _mm256_unpacklo_epi32(t, r23);
	r23 = _mm256_unpackhi_epi32(t, r23);
	t = r27;
	r27 = _mm256_unpacklo_epi32(t, r31);
	r31 = _mm256_unpackhi_epi32(t, r31);

	t = r0;
	r0 = _mm256_unpacklo_epi64(t, r8);
	r8 = _mm256_unpackhi_epi64(t, r8);
	t = r16;
	r16 = _mm256_unpacklo_epi64(t, r24);
	r24 = _mm256_unpackhi_epi64(t, r24);
	t = r4;
	r4 = _mm256_unpacklo_epi64(t, r12);
	r12 = _mm256_unpackhi_epi64(t, r12);
	t = r20;
	r20 = _mm256_unpacklo_epi64(t, r28);
	r28 = _mm256_unpackhi_epi64(t, r28);
	t = r2;
	r2 = _mm256_unpacklo_epi64(t, r10);
	r10 = _mm256_unpackhi_epi64(t, r10);
	t = r18;
	r18 = _mm256_unpacklo_epi64(t, r26);
	r26 = _mm256_unpackhi_epi64(t, r26);
	t = r6;
	r6 = _mm256_unpacklo_epi64(t, r14);
	r14 = _mm256_unpackhi_epi64(t, r14);
	t = r22;
	r22 = _mm256_unpacklo_epi64(t, r30);
	r30 = _mm256_unpackhi_epi64(t, r30);
	t = r1;
	r1 = _mm256_unpacklo_epi64(t, r9);
	r9 = _mm256_unpackhi_epi64(t, r9);
	t = r17;
	r17 = _mm256_unpacklo_epi64(t, r25);
	r25 = _mm256_unpackhi_epi64(t, r25);
	t = r5;
	r5 = _mm256_unpacklo_epi64(t, r13);
	r13 = _mm256_unpackhi_epi64(t, r13);
	t = r21;
	r21 = _mm256_unpacklo_epi64(t, r29);
	r29 = _mm256_unpackhi_epi64(t, r29);
	t = r3;
	r3 = _mm256_unpacklo_epi64(t, r11);
	r11 = _mm256_unpackhi_epi64(t, r11);
	t = r19;
	r19 = _mm256_unpacklo_epi64(t, r27);
	r27 = _mm256_unpackhi_epi64(t, r27);
	t = r7;
	r7 = _mm256_unpacklo_epi64(t, r15);
	r15 = _mm256_unpackhi_epi64(t, r15);
	t = r23;
	r23 = _mm256_unpacklo_epi64(t, r31);
	r31 = _mm256_unpackhi_epi64(t, r31);

	t = r0;
	r0 = _mm256_permute2x128_si256(t, r16, 0x20);
	r16 = _mm256_permute2x128_si256(t, r16, 0x31);
	t = r8;
	r8 = _mm256_permute2x128_si256(t, r24, 0x20);
	r24 = _mm256_permute2x128_si256(t, r24, 0x31);
	t = r4;
	r4 = _mm256_permute2x128_si256(t, r20, 0x20);
	r20 = _mm256_permute2x128_si256(t, r20, 0x31);
	t = r12;
	r12 = _mm256_permute2x128_si256(t, r28, 0x20);
	r28 = _mm256_permute2x128_si256(t, r28, 0x31);
	t = r2;
	r2 = _mm256_permute2x128_si256(t, r18, 0x20);
	r18 = _mm256_permute2x128_si256(t, r18, 0x31);
	t = r10;
	r10 = _mm256_permute2x128_si256(t, r26, 0x20);
	r26 = _mm256_permute2x128_si256(t, r26, 0x31);
	t = r6;
	r6 = _mm256_permute2x128_si256(t, r22, 0x20);
	r22 = _mm256_permute2x128_si256(t, r22, 0x31);
	t = r14;
	r14 = _mm256_permute2x128_si256(t, r30, 0x20);
	r30 = _mm256_permute2x128_si256(t, r30, 0x31);
	t = r1;
	r1 = _mm256_permute2x128_si256(t, r17, 0x20);
	r17 = _mm256_permute2x128_si256(t, r17, 0x31);
	t = r9;
	r9 = _mm256_permute2x128_si256(t, r25, 0x20);
	r25 = _mm256_permute2x128_si256(t, r25, 0x31);
	t = r5;
	r5 = _mm256_permute2x128_si256(t, r21, 0x20);
	r21 = _mm256_permute2x128_si256(t, r21, 0x31);
	t = r13;
	r13 = _mm256_permute2x128_si256(t, r29, 0x20);
	r29 = _mm256_permute2x128_si256(t, r29, 0x31);
	t = r3;
	r3 = _mm256_permute2x128_si256(t, r19, 0x20);
	r19 = _mm256_permute2x128_si256(t, r19, 0x31);
	t = r11;
	r11 = _mm256_permute2x128_si256(t, r27, 0x20);
	r27 = _mm256_permute2x128_si256(t, r27, 0x31);
	t = r7;
	r7 = _mm256_permute2x128_si256(t, r23, 0x20);
	r23 = _mm256_permute2x128_si256(t, r23, 0x31);
	t = r15;
	r15 = _mm256_permute2x128_si256(t, r31, 0x20);
	r31 = _mm256_permute2x128_si256(t, r31, 0x31);
	
	__m256i* ptr = (__m256i*)out;
	_mm256_store_si256(ptr++, r0);
	_mm256_store_si256(ptr++, r16);
	_mm256_store_si256(ptr++, r8);
	_mm256_store_si256(ptr++, r24);
	_mm256_store_si256(ptr++, r4);
	_mm256_store_si256(ptr++, r20);
	_mm256_store_si256(ptr++, r12);
	_mm256_store_si256(ptr++, r28);
	_mm256_store_si256(ptr++, r2);
	_mm256_store_si256(ptr++, r18);
	_mm256_store_si256(ptr++, r10);
	_mm256_store_si256(ptr++, r26);
	_mm256_store_si256(ptr++, r6);
	_mm256_store_si256(ptr++, r22);
	_mm256_store_si256(ptr++, r14);
	_mm256_store_si256(ptr++, r30);
	_mm256_store_si256(ptr++, r1);
	_mm256_store_si256(ptr++, r17);
	_mm256_store_si256(ptr++, r9);
	_mm256_store_si256(ptr++, r25);
	_mm256_store_si256(ptr++, r5);
	_mm256_store_si256(ptr++, r21);
	_mm256_store_si256(ptr++, r13);
	_mm256_store_si256(ptr++, r29);
	_mm256_store_si256(ptr++, r3);
	_mm256_store_si256(ptr++, r19);
	_mm256_store_si256(ptr++, r11);
	_mm256_store_si256(ptr++, r27);
	_mm256_store_si256(ptr++, r7);
	_mm256_store_si256(ptr++, r23);
	_mm256_store_si256(ptr++, r15);
	_mm256_store_si256(ptr, r31);
}

#endif
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
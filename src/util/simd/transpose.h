/****
DIAMOND protein aligner
Copyright (C) 2013-2019 Benjamin Buchfink <buchfink@gmail.com>

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

#ifndef TRANSPOSE_H_
#define TRANSPOSE_H_

#include "../simd.h"

inline void merge(__m128i &first, __m128i &second, __m128i &temp)
{
	temp = first;
	first = _mm_unpacklo_epi8(first, second);
	temp = _mm_unpackhi_epi8(temp, second);
	second = temp;
}

inline void merge2(__m128i &col0, __m128i &col1, __m128i &col2, __m128i &col3, __m128i &temp)
{
	merge(col0, col2, temp);
	merge(col1, col3, temp);
	merge(col0, col1, temp);
	merge(col2, col3, temp);
}

inline void merge3(__m128i &col0, __m128i &col1, __m128i &col2, __m128i &col3, __m128i &col4, __m128i &col5, __m128i &col6, __m128i &col7, __m128i &tmp)
{
	merge2(col0, col2, col4, col6, tmp);
	merge2(col1, col3, col5, col7, tmp);
	merge(col0, col1, tmp);
	merge(col2, col3, tmp);
	merge(col4, col5, tmp);
	merge(col6, col7, tmp);
}

inline void merge_fast(__m128i &in0, __m128i &in1, __m128i &tmp)
{
	tmp = in0;
	in0 = _mm_unpacklo_epi8(in0, in1);
	tmp = _mm_unpackhi_epi8(tmp, in1);
}

inline void merge_fast_write(__m128i &reg1, __m128i &reg2, int off1, int off2, __m128i &tmp, __m128i *outreg)
{
	merge_fast(reg1, reg2, tmp);
	_mm_store_si128(outreg+off1, reg1);
	_mm_store_si128(outreg+off2, tmp);
}

inline void merge2_fast(__m128i &in0, __m128i &in1, __m128i &in2, __m128i &in3, __m128i &tmp)
{
	merge_fast(in0, in2, tmp);
	merge_fast(in1, in3, in2);
	merge_fast(in0, in1, in3);
	merge_fast(tmp, in2, in1);
}

inline void merge3_fast(__m128i &in0, __m128i &in1, __m128i &in2, __m128i &in3, __m128i &in4, __m128i &in5, __m128i &in6, __m128i &in7, __m128i &tmp)
{
	merge2_fast(in0, in2, in4, in6, tmp);
	merge2_fast(in1, in3, in5, in7, in4);
	merge_fast(in0, in1, in5);
	merge_fast(in6, in7, in1);
	merge_fast(tmp, in4, in7);
	merge_fast(in2, in3, in4);
}

inline void transpose(char* in, char* out, int hIn) {
	__m128i xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8, xmm9, xmm10, xmm11, xmm12, xmm13, xmm14, xmm15;
	xmm0 = _mm_load_si128((__m128i*)in); in += 16;
	xmm1 = _mm_load_si128((__m128i*)in); in += 16;
	xmm2 = _mm_load_si128((__m128i*)in); in += 16;
	xmm3 = _mm_load_si128((__m128i*)in); in += 16;
	xmm4 = _mm_load_si128((__m128i*)in); in += 16;
	xmm5 = _mm_load_si128((__m128i*)in); in += 16;
	xmm6 = _mm_load_si128((__m128i*)in); in += 16;
	xmm7 = _mm_load_si128((__m128i*)in); in += 16;
	xmm8 = _mm_load_si128((__m128i*)in); in += 16;
	xmm9 = _mm_load_si128((__m128i*)in); in += 16;
	xmm10 = _mm_load_si128((__m128i*)in); in += 16;
	xmm11 = _mm_load_si128((__m128i*)in); in += 16;
	xmm12 = _mm_load_si128((__m128i*)in); in += 16;
	xmm13 = _mm_load_si128((__m128i*)in); in += 16;
	xmm14 = _mm_load_si128((__m128i*)in); in += 16;

	merge3_fast(xmm0, xmm2, xmm4, xmm6, xmm8, xmm10, xmm12, xmm14, xmm15);
	xmm6 = xmm15;

	_mm_store_si128((__m128i*)out, xmm0);
	xmm15 = _mm_load_si128((__m128i*)in);

	merge3_fast(xmm1, xmm3, xmm5, xmm7, xmm9, xmm11, xmm13, xmm15, xmm0);
	xmm7 = xmm0;

	merge_fast_write(xmm8, xmm9, 14, 15, xmm0, (__m128i*)out);
	xmm0 = _mm_load_si128((__m128i*)out);

	merge_fast_write(xmm4, xmm5, 12, 13, xmm9, (__m128i*)out);
	merge_fast_write(xmm14, xmm15, 10, 11, xmm9, (__m128i*)out);
	merge_fast_write(xmm6, xmm7, 8, 9, xmm9, (__m128i*)out);
	merge_fast_write(xmm2, xmm3, 6, 7, xmm9, (__m128i*)out);
	merge_fast_write(xmm12, xmm13, 4, 5, xmm9, (__m128i*)out);
	merge_fast_write(xmm10, xmm11, 2, 3, xmm9, (__m128i*)out);
	merge_fast_write(xmm0, xmm1, 0, 1, xmm9, (__m128i*)out);
}

#endif /* TRANSPOSE_H_ */

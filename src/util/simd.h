/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#ifndef SIMD_H_
#define SIMD_H_

#include "system.h"

#ifdef _MSC_VER
#define __MMX__
#define __SSE__
#define __SSE2__
#define __SSSE3__
#define __SSE4_1__
#endif

#ifdef __SSSE3__
#include <tmmintrin.h>
#endif
#ifdef __SSE4_1__
#include <smmintrin.h>
#endif
#ifdef __MMX__
#include <mmintrin.h>
#endif
#ifdef __SSE__
#include <xmmintrin.h>
#endif
#ifdef __SSE2__
#include <emmintrin.h>
#endif

inline void check_simd()
{
#ifdef __SSSE3__
	int info[4];
	cpuid(info, 0);
	int nids = info[0];
	if (nids >= 1) {
		cpuid(info, 1);
		if ((info[2] & (1 << 9)) == 0)
			throw std::runtime_error("CPU does not support SSSE3. Please try to run the diamond-sse2 file contained in the binary package or compile the software from source.");
	}
	else
		throw std::runtime_error("Incompatible CPU type. Please try to run the diamond-sse2 file contained in the binary package or compile the software from source.");
#endif
}

#endif
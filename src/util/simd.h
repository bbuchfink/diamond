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

#ifndef SIMD_H_
#define SIMD_H_

#include "system.h"

#ifdef _MSC_VER
#define __SSE2__
#define __SSSE3__
#define __SSE4_1__
#define __POPCNT__
#endif

#ifdef __SSSE3__
#include <tmmintrin.h>
#endif
#ifdef __SSE2__
#include <emmintrin.h>
#endif
#ifdef __SSE4_1__
#include <smmintrin.h>
#endif

namespace SIMD {

enum class Arch { Generic, SSE4_2 };
enum Flags { SSSE3 = 1, POPCNT = 2, SSE4_2 = 4 };
extern Arch arch;

#define DISPATCH(x) switch(SIMD::arch) {\
case Arch::Generic: DISPATCH_GENERIC::x; break;\
case Arch::SSE4_2: DISPATCH_SSE4_2::x;}

#define DECL_DISPATCH(x) namespace ARCH_GENERIC { x; }\
namespace ARCH_SSE4_2 { x; }

void init();

};

void check_simd();
void simd_messages();

#endif
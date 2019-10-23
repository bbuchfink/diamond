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

#include <string>
#include "system.h"

#ifdef _MSC_VER
#define __SSE2__
#ifdef DISPATCH_ARCH
#if DISPATCH_ARCH == ARCH_SSE4_1
#define __SSSE3__
#define __SSE4_1__
#define __POPCNT__
#endif
#endif
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
#ifdef _MSC_VER
#include <intrin.h>
#endif

namespace SIMD {

enum class Arch { Generic, SSE4_1 };
enum Flags { SSSE3 = 1, POPCNT = 2, SSE4_1 = 4 };
extern Arch arch;
extern int flags;

#define DISPATCH(ns, f) (SIMD::arch == SIMD::Arch::SSE4_1 ? ns ## ARCH_SSE4_1::f : ns ## ARCH_GENERIC::f)

#define DECL_DISPATCH(f) namespace ARCH_GENERIC { f; }\
namespace ARCH_SSE4_1 { f; }

void init();
std::string features();

};

#endif
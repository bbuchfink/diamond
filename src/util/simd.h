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
#include <functional>
#include "system.h"

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

enum class Arch { None, Generic, SSE4_1 };
enum Flags { SSSE3 = 1, POPCNT = 2, SSE4_1 = 4 };
Arch arch();

#ifdef __SSE__
#define DECL_DISPATCH(ret, name, param) namespace ARCH_GENERIC { ret name param; }\
namespace ARCH_SSE4_1 { ret name param; }\
inline std::function<decltype(ARCH_GENERIC::name)> dispatch_target_##name() {\
switch(SIMD::arch()) {\
case SIMD::Arch::SSE4_1: return ARCH_SSE4_1::name;\
default: return ARCH_GENERIC::name;\
}}\
const std::function<decltype(ARCH_GENERIC::name)> name = dispatch_target_##name();
#else
#define DECL_DISPATCH(ret, name, param) namespace ARCH_GENERIC { ret name param; }\
namespace ARCH_SSE4_1 { ret name param; }\
inline std::function<decltype(ARCH_GENERIC::name)> dispatch_target_##name() {\
return ARCH_GENERIC::name;\
}\
const std::function<decltype(ARCH_GENERIC::name)> name = dispatch_target_##name();
#endif

std::string features();

};

#endif

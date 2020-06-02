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

#include <vector>
#include "simd.h"
#include "util.h"

#ifdef _WIN32
#define cpuid(info,x)    __cpuidex(info,x,0)
#else
#ifdef __SSE2__
inline void cpuid(int CPUInfo[4], int InfoType) {
	__asm__ __volatile__(
		"cpuid":
	"=a" (CPUInfo[0]),
		"=b" (CPUInfo[1]),
		"=c" (CPUInfo[2]),
		"=d" (CPUInfo[3]) :
		"a" (InfoType), "c" (0)
		);
}
#endif
#endif

namespace SIMD {

int flags = 0;

Arch init_arch() {
#ifdef __SSE2__
	int info[4];
	cpuid(info, 0);
	int nids = info[0];
	if (nids >= 1) {
		cpuid(info, 1);
	}
	else
		throw std::runtime_error("Incompatible CPU type. Please try to compile the software from source.");

	if ((info[2] & (1 << 9)) != 0)
		flags |= SSSE3;
	if ((info[2] & (1 << 23)) != 0)
		flags |= POPCNT;
	if ((info[2] & (1 << 19)) != 0)
		flags |= SSE4_1;
	if (nids >= 7) {
		cpuid(info, 7);
		if ((info[1] & (1 << 5)) != 0)
			flags |= AVX2;
	}
#endif

#ifdef __SSSE3__
	if ((flags & SSSE3) == 0)
		throw std::runtime_error("CPU does not support SSSE3. Please compile the software from source.");
#endif
#ifdef __POPCNT__
	if ((flags & POPCNT) == 0)
		throw std::runtime_error("CPU does not support POPCNT. Please compile the software from source.");
#endif
#ifdef __SSE4_1__
	if ((flags & SSE4_1) == 0)
		throw std::runtime_error("CPU does not support SSE4.1. Please compile the software from source.");
#endif
#ifdef __AVX2__
	if ((flags & AVX2) == 0)
		throw std::runtime_error("CPU does not support AVX2. Please compile the software from source.");
#endif

	if ((flags & SSSE3) && (flags & POPCNT) && (flags & SSE4_1) && (flags & AVX2))
		return Arch::AVX2;
	if ((flags & SSSE3) && (flags & POPCNT) && (flags & SSE4_1))
		return Arch::SSE4_1;
	else
		return Arch::Generic;
}

Arch arch() {
	static Arch a = Arch::None;
	return a == Arch::None ? (a = init_arch()) : a;
}

std::string features() {
	std::vector<std::string> r;
	if (flags & SSSE3)
		r.push_back("ssse3");
	if (flags & POPCNT)
		r.push_back("popcnt");
	if (flags & SSE4_1)
		r.push_back("sse4.1");
	if (flags & AVX2)
		r.push_back("avx2");
	return r.empty() ? "None" : join(" ", r);
}

}
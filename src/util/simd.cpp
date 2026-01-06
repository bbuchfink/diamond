/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

#include <stdexcept>
#include <vector>
#include "simd.h"
#include "util/string/string.h"

#ifdef HAVE_GETAUXVAL
#include <sys/auxv.h>
#endif

using std::string;

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
#if defined(__aarch64__) | defined(__arm__)

#if defined(__aarch64__)
	flags |= NEON;
#elif defined(HAVE_MFPU_NEON) & defined(HAVE_GETAUXVAL)
	unsigned long hwcaps = getauxval(AT_HWCAP);
	if ((hwcaps & HWCAP_ARM_NEON) != 0)
		flags |= NEON;
#endif

	if ((flags & NEON) != 0)
		return Arch::NEON;

#else

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
#ifdef WITH_AVX512
		if ((info[1] & (1 << 16)) != 0 && (info[1] & (1 << 30)) != 0)
			flags |= AVX512;
#endif
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

	if (flags & AVX512)
		return Arch::AVX512;
	if ((flags & SSSE3) && (flags & POPCNT) && (flags & SSE4_1) && (flags & AVX2))
		return Arch::AVX2;
	if ((flags & SSSE3) && (flags & POPCNT) && (flags & SSE4_1))
		return Arch::SSE4_1;
#endif
	return Arch::Generic;
}

Arch arch() {
	static Arch a = Arch::None;
	return a == Arch::None ? (a = init_arch()) : a;
}

string features() {
	init_arch();
	std::vector<string> r;
	if (flags & NEON)
		r.push_back("neon");
	if (flags & SSSE3)
		r.push_back("ssse3");
	if (flags & POPCNT)
		r.push_back("popcnt");
	if (flags & SSE4_1)
		r.push_back("sse4.1");
	if (flags & AVX2)
		r.push_back("avx2");
	return r.empty() ? "None" : join(" ", r.begin(), r.end());
}

}
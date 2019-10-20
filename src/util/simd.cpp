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

#include "simd.h"
#include "log_stream.h"

#ifdef _WIN32
#define cpuid(info,x)    __cpuidex(info,x,0)
#else
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

void check_simd()
{
#ifdef __SSE2__
	int info[4];
	cpuid(info, 0);
	int nids = info[0];
	if (nids >= 1) {
		cpuid(info, 1);
	}
	else
		throw std::runtime_error("Incompatible CPU type. Please try to compile the software from source.");
#endif
#ifdef __SSSE3__
	if ((info[2] & (1 << 9)) == 0)
		throw std::runtime_error("CPU does not support SSSE3. Please compile the software from source.");
#endif
#ifdef __POPCNT__
	if ((info[2] & (1 << 23)) == 0)
		throw std::runtime_error("CPU does not support POPCNT. Please compile the software from source.");
#endif
#ifdef __SSE4_1__
	if ((info[2] & (1 << 19)) == 0)
		throw std::runtime_error("CPU does not support SSE4.1. Please compile the software from source.");
#endif
}

void simd_messages()
{
#ifdef __SSSE3__
	verbose_stream << "SSSE3 enabled." << endl;
#endif
#ifdef __SSE4_1__
	verbose_stream << "SSE4.1 enabled." << endl;
#endif
#ifdef __POPCNT__
	verbose_stream << "POPCNT enabled." << endl;
#endif
}

namespace SIMD {

Arch arch;

void init() {
#ifdef __SSE2__
	int info[4];
	cpuid(info, 0);
	int nids = info[0];
	if (nids >= 1) {
		cpuid(info, 1);
	}
	else
		throw std::runtime_error("Incompatible CPU type. Please try to compile the software from source.");
#endif
}

}
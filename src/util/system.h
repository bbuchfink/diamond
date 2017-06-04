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

#ifndef SYSTEM_H_
#define SYSTEM_H_

#include <stdexcept>

#ifdef _WIN32
#define cpuid(info,x)    __cpuidex(info,x,0)
#else
inline void cpuid(int CPUInfo[4],int InfoType) {
    __asm__ __volatile__ (
        "cpuid":
        "=a" (CPUInfo[0]),
        "=b" (CPUInfo[1]),
        "=c" (CPUInfo[2]),
        "=d" (CPUInfo[3]) :
        "a" (InfoType), "c" (0)
    );
}
#endif

#ifdef _MSC_VER

#define PACKED_ATTRIBUTE

#else

#define PACKED_ATTRIBUTE __attribute__((packed))
#define FTELL(x) ftell(x)
#define FSEEK(x,y,z) fseek(x,y,z)

#endif

#ifndef _MSC_VER

#include <sys/resource.h>

inline void set_max_open_files(unsigned n)
{
	rlimit rlp;
	if (getrlimit(RLIMIT_NOFILE, &rlp) != 0)
		throw std::runtime_error("Error executing getrlimit.");
	if (rlp.rlim_max < n)
		throw std::runtime_error("Open files hard limit is too low. Set lower value for --bin parameter.");
	if (rlp.rlim_cur < n) {
		rlp.rlim_cur = n;
		if (setrlimit(RLIMIT_NOFILE, &rlp) != 0)
			throw std::runtime_error("Error executing setrlimit.");
	}
	//std::cout << "Soft limit = " << rlp.rlim_cur << " Hard limit = " << rlp.rlim_max << std::endl;
}

#else

inline void set_max_open_files(unsigned n)
{}

#endif

#endif /* SYSTEM_H_ */

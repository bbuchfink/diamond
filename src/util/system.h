/****
Copyright (c) 2014, University of Tuebingen
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****
Author: Benjamin Buchfink
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
#define FTELL(x) _ftelli64(x)
#define FSEEK(x,y,z) _fseeki64(x,y,z)

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

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

#pragma once
#include <cstdio>
#include <exception>
#include <string>

#ifdef _MSC_VER

typedef __int64 ssize_t;
#define STDIN_FILENO _fileno(stdin)
#define STDOUT_FILENO _fileno(stdout)
#define PACKED_ATTRIBUTE
#define FORCE_INLINE __forceinline static
#define FLATTEN
#define UNLINK _unlink
#define NO_INLINE

#else

#ifdef __sparc__
#define PACKED_ATTRIBUTE
#else
#define PACKED_ATTRIBUTE __attribute__((packed))
#endif
#define FORCE_INLINE __attribute__((always_inline)) static inline
#define FLATTEN __attribute__((flatten))
#define UNLINK unlink
#define NO_INLINE __attribute__ ((noinline))

#endif

#ifdef __linux__
#define POSIX_OPEN(x,y,z) open64(x,y,z)
#define POSIX_OPEN2(x,y) open64(x,y)
#else
#define POSIX_OPEN(x,y,z) open(x,y,z)
#define POSIX_OPEN2(x,y) open(x,y)
#endif
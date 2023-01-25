/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen
						
Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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

#pragma once
#include <stdexcept>

#ifdef _MSC_VER

typedef __int64 ssize_t;
#define STDIN_FILENO _fileno(stdin)
#define STDOUT_FILENO _fileno(stdout)
#define PACKED_ATTRIBUTE
#define FORCE_INLINE __forceinline static
#define FLATTEN
#define UNLINK _unlink

#else

#define PACKED_ATTRIBUTE __attribute__((packed))
#define FORCE_INLINE __attribute__((always_inline)) static inline
#define FLATTEN __attribute__((flatten))
#define UNLINK unlink

#endif

#ifdef __linux__
#define POSIX_OPEN(x,y,z) open64(x,y,z)
#define POSIX_OPEN2(x,y) open64(x,y)
#else
#define POSIX_OPEN(x,y,z) open(x,y,z)
#define POSIX_OPEN2(x,y) open(x,y)
#endif
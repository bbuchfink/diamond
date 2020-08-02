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

#define PACKED_ATTRIBUTE

#else

#define PACKED_ATTRIBUTE __attribute__((packed))

#endif

#ifdef __linux__
#define POSIX_OPEN(x,y,z) open64(x,y,z)
#define POSIX_OPEN2(x,y) open64(x,y)
#else
#define POSIX_OPEN(x,y,z) open(x,y,z)
#define POSIX_OPEN2(x,y) open(x,y)
#endif

#if defined(__s390x__) && defined(__clang__)
#define TLS_FIX_S390X
#define TLS_FIX_S390X_MOVE(x) x
#else
#define TLS_FIX_S390X thread_local
#define TLS_FIX_S390X_MOVE(x) std::move(x)
#endif
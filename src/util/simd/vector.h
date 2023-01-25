/****
DIAMOND protein aligner
Copyright (C) 2020 Max Planck Society for the Advancement of Science e.V.

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

#ifndef SIMD_VECTOR_H_
#define SIMD_VECTOR_H_

#include "../simd.h"

#if ARCH_ID == 3
#include "vector8_avx512.h"
#elif ARCH_ID == 2
#include "vector8_avx2.h"
#elif defined(__SSE2__)
#include "vector8_sse.h"
#else
#include "vector_generic.h"
#endif

#endif
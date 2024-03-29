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

#pragma once
#include "../simd.h"

#if defined(__SSE2__) | defined(__ARM_NEON)
#include "transpose16x16.h"
#endif

#if ARCH_ID == 2
#include "transpose32x32.h"
#endif

#if ARCH_ID == 3
static inline void transpose(const signed char** data, size_t n, signed char* out, const __m512i&) {}
#endif

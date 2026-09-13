/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once
#include "../simd.h"

#if defined(__SSE2__) | defined(__ARM_NEON)
#include "transpose16x16.h"
#endif

#if ARCH_AVX2_KERNELS
#include "transpose32x32.h"
#endif

#if !defined(__SSE2__) && !defined(__ARM_NEON)

#include <stddef.h>
#include <stdint.h>

static inline void transpose_offset(const int16_t** data, size_t n, ptrdiff_t offset, int16_t* out, int16_t) {
	*out = (*data)[offset];
}

#endif
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
#include "util/enum.h"

namespace DP {

	enum class Flags { NONE = 0, PARALLEL = 1, FULL_MATRIX = 2, SEMI_GLOBAL = 4 };

	DEFINE_ENUM_FLAG_OPERATORS(Flags)

}

enum class HspValues : unsigned {
	NONE = 0,
	TRANSCRIPT = 1,
	QUERY_START = 1 << 1,
	QUERY_END = 1 << 2,
	TARGET_START = 1 << 3,
	TARGET_END = 1 << 4,
	IDENT = 1 << 5,
	LENGTH = 1 << 6,
	MISMATCHES = 1 << 7,
	GAP_OPENINGS = 1 << 8,
	GAPS = IDENT | LENGTH | MISMATCHES,
	QUERY_COORDS = QUERY_START | QUERY_END,
	TARGET_COORDS = TARGET_START | TARGET_END,
	COORDS = QUERY_COORDS | TARGET_COORDS
};

DEFINE_ENUM_FLAG_OPERATORS(HspValues)

static inline bool have_coords(const HspValues v) {
	return flag_any(v, HspValues::TRANSCRIPT) || flag_all(v, HspValues::QUERY_COORDS | HspValues::TARGET_COORDS);
}
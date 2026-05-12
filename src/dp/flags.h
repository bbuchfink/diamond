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
#include "util/enum.h"

namespace DP {

	enum class Flags { NONE = 0, PARALLEL = 1, FULL_MATRIX = 2, SEMI_GLOBAL = 4 };

	DEFINE_ENUM_OPERATORS(Flags)

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

DEFINE_ENUM_OPERATORS(HspValues)

static inline bool have_coords(const HspValues v) {
	return flag_any(v, HspValues::TRANSCRIPT) || flag_all(v, HspValues::QUERY_COORDS | HspValues::TARGET_COORDS);
}
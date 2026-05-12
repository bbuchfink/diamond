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

namespace Output {

enum class Flags : int {
	NONE = 0,
	FULL_TITLES = 1,
	ALL_SEQIDS = 1 << 1,
	TARGET_SEQS = 1 << 2,
	SELF_ALN_SCORES = 1 << 3,
    IS_STRING = 1 << 4,
    IS_ARRAY = 1 << 5,
	SSEQID = 1 << 6,
	DEFAULT_REPORT_UNALIGNED = 1 << 7,
	NO_REALIGN = 1 << 8
};

DEFINE_ENUM_OPERATORS(Flags)

}
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
#include <stdint.h>
#include <stddef.h>
#include <limits>
#include "basic/value.h"

namespace Cluster { namespace Multinode {

constexpr uint64_t LEN_SORT_BLOCK_SEQ_LIMIT = (uint64_t)std::numeric_limits<BlockId>::max() - 1;
constexpr uint64_t LEN_SORT_BLOCK_RAW_LIMIT = (uint64_t(1) << 40) - 1;
constexpr uint64_t LEN_SORT_BLOCK_RAW_PADDING = 256;

bool can_add_to_len_sorted_block(
	uint64_t block_letters,
	uint64_t block_seqs,
	uint64_t seq_len,
	uint64_t block_letter_limit,
	uint64_t block_seq_limit = LEN_SORT_BLOCK_SEQ_LIMIT,
	uint64_t block_raw_limit = LEN_SORT_BLOCK_RAW_LIMIT);

double block_combo_chunk_size(size_t db_file_size, size_t query_file_size);

} }

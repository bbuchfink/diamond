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
#include <vector>
#include "../basic/packed_loc.h"
#include "../masking/def.h"
#include "basic/sequence.h"
#include "util/data_structures/bit_vector.h"

enum class SeedEncoding { SPACED_FACTOR, HASHED, CONTIGUOUS };

struct NoFilter
{
	bool contains(uint64_t seed, uint64_t shape) const
	{
		return true;
	}
};

extern NoFilter no_filter;

#ifndef __sparc__
#pragma pack(1)
#endif
struct PackedLocId {
	PackedLocId() {}
	PackedLocId(PackedLoc pos) :
		pos(pos) {}
	PackedLocId(PackedLoc pos, uint32_t block_id) :
		pos(pos),
		block_id(block_id)
	{}
	operator uint64_t() const {
		return (uint64_t)pos;
	}
	PackedLoc pos;
	uint32_t block_id;
} PACKED_ATTRIBUTE;
#pragma pack()

static inline uint32_t block_id(PackedLocId i) {
	return i.block_id;
}

static inline uint32_t block_id(PackedLoc i) {
	throw std::runtime_error("Unsupported");
}

struct EnumCfg {
	const std::vector<uint32_t>* partition;
	int shape_begin, shape_end;
	const SeedEncoding code;
	const std::vector<bool>* const skip;
	const bool filter_masked_seeds, mask_seeds;
	const double seed_cut;
	const MaskingAlgo soft_masking;
	const Loc minimizer_window;
	const bool filter_low_complexity_seeds, mask_low_complexity_seeds;
	const Loc sketch_size;
	const std::vector<BitVector>* const skip_seed_positions;
};

struct SeqInfo {
	BlockId block_id;
	OId oid;
	const char* title, * qual;
	Loc len;
	Sequence source_seq, mate_seq;
};
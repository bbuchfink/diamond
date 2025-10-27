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
#include <stdint.h>
#include <vector>
#include "../basic/packed_loc.h"
#include "../masking/def.h"
#include "basic/sequence.h"

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
};

struct SeqInfo {
	BlockId block_id;
	OId oid;
	const char* title, * qual;
	Loc len;
	Sequence source_seq, mate_seq;
};
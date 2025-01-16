#pragma once
#include <stdint.h>
#include <vector>
#include "../basic/packed_loc.h"
#include "../masking/def.h"

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
	const std::vector<BlockId>* partition;
	size_t shape_begin, shape_end;
	const SeedEncoding code;
	const std::vector<bool>* const skip;
	const bool filter_masked_seeds, mask_seeds;
	const double seed_cut;
	const MaskingAlgo soft_masking;
	const Loc minimizer_window;
	const bool filter_low_complexity_seeds, mask_low_complexity_seeds;
	const Loc sketch_size;
};
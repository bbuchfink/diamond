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

#ifdef KEEP_TARGET_ID
struct SeedLoc {
	SeedLoc() {}
	SeedLoc(PackedLoc pos) :
		pos(pos) {}
	SeedLoc(PackedLoc pos, uint32_t block_id) :
		pos(pos),
		block_id(block_id)
	{}
	operator uint64_t() const {
		return (uint64_t)pos;
	}
	PackedLoc pos;
	uint32_t block_id;
};
#else
using SeedLoc = PackedLoc;
#endif

struct EnumCfg {
	const std::vector<size_t>* partition;
	const size_t shape_begin, shape_end;
	const SeedEncoding code;
	const std::vector<bool>* const skip;
	const bool filter_masked_seeds, mask_seeds;
	const double seed_cut;
	const MaskingAlgo soft_masking;
	const Loc minimizer_window;
};

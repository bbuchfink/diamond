#pragma once
#include "../util/enum.h"

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
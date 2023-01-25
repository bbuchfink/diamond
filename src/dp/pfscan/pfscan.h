#pragma once

#include "../basic/statistics.h"
#include "../../basic/match.h"
#include "../score_profile.h"

namespace DP { namespace PrefixScan {

struct Config {
	Sequence query, target;
	const char* query_seqid, *target_seqid;
	Loc d_begin, d_end;
	Interval hint_target_range;
	const int16_t* const* query_profile, * const* query_profile_rev;
	const int8_t* const* query_profile8, * const* query_profile_rev8;
	Statistics& stats;
	Score score_bias, xdrop, score_hint;

	void adjust_band(Loc channels) {
		const Loc b = std::max((d_end - d_begin + channels - 1) / channels * channels, channels);
		d_begin -= (b - d_end + d_begin) / 2;
		d_end = d_begin + b;
	}
	Loc band() const {
		return d_end - d_begin;
	}
};

DECL_DISPATCH(Hsp, align16, (const Config& cfg))
DECL_DISPATCH(Hsp, align8, (const Config& cfg))
DECL_DISPATCH(Hsp, align_anchored, (const Anchor& anchor, const Config& cfg))

}}
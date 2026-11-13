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

#include <algorithm>
#include <cassert>
#include <vector>
#include "dp/dp.h"
#include "align/def.h"
#include "dp/anchored_swipe/anchored.h"
#include "dp/long_profile/score_profile.h"
#include "util/simd/dispatch.h"
#include "util/geo/geo.h"
#include "util/log_stream.h"
#include "util/memory/mem_profile.h"
#include "stats/score_matrix.h"
#include "anchored_swipe.h"

namespace MultipassPipeline { namespace DISPATCH_ARCH {

static_assert(DP::AnchoredSwipe::MATRIX_ROW_SCORES, "The kernel is only supported in MATRIX_ROW_SCORES mode.");

using Target = DP::AnchoredSwipe::Target<int16_t>;
using ScoreVector = ::DISPATCH_ARCH::ScoreVector<int16_t, 0>;

// Returns the DP target of an extension, without its query letters, which are set once it is loaded into a
// channel of the kernel.
static Target make_target(const Extension& e, int64_t target_idx, Loc band_max) {
	const ExtensionMatrix m = e.matrix();
	assert(m.target_len > 0);
	assert(m.d_end >= 1 && m.band() >= 1 && m.band() <= band_max);
	const Sequence& target = e.pair->target;
	const bool reverse = e.direction == LEFT;
	// The target letters of the sub-matrix run from the first cell of the extension into its direction.
	const Sequence seq = reverse ? target.subseq(e.j + 1 - m.target_len, e.j + 1) : target.subseq(e.j, e.j + m.target_len);
	Target t(seq, m.d_begin, m.d_end, m.query_start, m.query_len, target_idx, reverse);
	t.profile = t.profile_rev = nullptr;
	t.query = t.query_rev = nullptr;
	return t;
}

// Writes seq framed by `padding` padding letters on each side to out, reversed if reverse is set, and returns the pointer to its first letter.
static const Letter* pad_query(Sequence seq, bool reverse, int64_t padding, Letter* out) {
	std::fill(out, out + padding, DP::AnchoredSwipe::PADDING_LETTER);
	Letter* begin = out + padding;
	if (reverse)
		std::reverse_copy(seq.data(), seq.end(), begin);
	else
		std::copy(seq.data(), seq.end(), begin);
	std::fill(begin + seq.length(), begin + seq.length() + padding, DP::AnchoredSwipe::PADDING_LETTER);
	return begin;
}

/* Target source of the kernel (see TargetArray in dp/anchored_swipe/anchored.h) feeding it a fixed range of extensions
   in the given order. Once they are used up, the kernel completes the channels still loaded and returns.
   The padded query of an extension is built when it is loaded into a channel and lives in the storage of that
   channel until the channel is loaded again. An extension whose pair carries a score table (a compositionally
   adjusted matrix) is computed with it, the others with the unadjusted matrix. A compact score table is expanded
   into a table of the channel, which lives as long as the padded query.
   Every extension is passed to the callback as soon as its channel is complete. */
class ExtensionSource {
public:

	ExtensionSource(const Extension* const* begin, const Extension* const* end, Loc band_max, const int8_t* common_scores, Statistics& stats, const ExtensionCallback& callback) :
		begin_(begin),
		pos_(begin),
		end_(end),
		band_max_(band_max),
		common_scores_(common_scores),
		stats_(stats),
		callback_(callback),
		padding_(std::max(DP::AnchoredSwipe::DISPATCH_ARCH::profile_padding<ScoreVector>(band_max), (int64_t)LongScoreProfile<int16_t>::DEFAULT_PADDING)),
		load_time_(0)
	{
		DP::AnchoredSwipe::make_score_table(score_matrix.matrix8(), 32, score_table_);
	}

	~ExtensionSource() {
		stats_.inc(Statistics::TIME_PROFILE, load_time_ / 1000);
	}

	bool next(int channel, Target& t) {
		if (pos_ == end_)
			return false;
		if ((size_t)channel >= channels_.size())
			channels_.resize(channel + 1);
		Channel& c = channels_[channel];
		c.extension = *pos_;
		t = make_target(*c.extension, pos_ - begin_, band_max_);
		++pos_;
		load(c, t);
		return true;
	}

	void finish(int channel, int16_t score, Loc query_end, Loc target_end) {
		callback_(*channels_[channel].extension, score, query_end, target_end);
	}

	DP::AnchoredSwipe::Options options() const {
		return DP::AnchoredSwipe::Options{ nullptr, nullptr, score_table_ };
	}

private:

	struct Channel {
		const Extension* extension = nullptr;
		std::vector<Letter> query;              // padded query of the extension, reversed for LEFT extensions
		// Score table (32x32) that the compact score table of the pair is expanded into. Holds the common scores
		// once it has been used.
		std::vector<int8_t> score_table;
	};

	// Builds the padded query of the extension of the channel, for the direction it extends into, and sets its
	// score table.
	void load(Channel& c, Target& t) {
		MEM_SCOPE("dp/query-profile");
		TaskTimer timer;
		const TargetPair& pair = *c.extension->pair;
		const bool reverse = c.extension->direction == LEFT;
		c.query.resize((size_t)pair.query.length() + 2 * padding_);
		(reverse ? t.query_rev : t.query) = pad_query(pair.query, reverse, padding_, c.query.data());
		if (pair.score_table_compact) {
			if (c.score_table.empty())
				c.score_table.assign(common_scores_, common_scores_ + 32 * 32);
			reinterpret_cast<const CompactScoreTable*>(pair.score_table)->expand(c.score_table.data());
			t.score_table = c.score_table.data();
		}
		else
			t.score_table = pair.score_table;
		load_time_ += timer.nanoseconds();
	}

	const Extension* const* const begin_;
	const Extension* const* pos_;              // next extension to load
	const Extension* const* const end_;
	const Loc band_max_;
	const int8_t* const common_scores_;
	Statistics& stats_;
	const ExtensionCallback& callback_;
	const int64_t padding_;
	alignas(32) int8_t score_table_[32 * 32];  // unadjusted matrix
	std::vector<Channel> channels_;
	int64_t load_time_;
};

void anchored_swipe(const Extension* const* begin, const Extension* const* end, Loc band_max, const int8_t* common_scores, Statistics& stats, const ExtensionCallback& callback) {
	MEM_SCOPE("dp/anchored-swipe");
	TaskTimer total;

	if (begin == end)
		return;

	{
		ExtensionSource source(begin, end, band_max, common_scores, stats, callback);
		const DP::AnchoredSwipe::Options options = source.options();

		TaskTimer timer;
		const DP::AnchoredSwipe::Stats kernel_stats = DP::AnchoredSwipe::DISPATCH_ARCH::smith_waterman<ScoreVector>(source, band_max, options);
#ifdef ANCHORED_SWIPE_STAT
		stats.inc(Statistics::ANCHORED_LANE_CELLS, kernel_stats.lane_cells);
		stats.inc(Statistics::ANCHORED_LANE_CELLS_IDLE, kernel_stats.idle_cells);
		stats.inc(Statistics::ANCHORED_LANE_CELLS_TARGET_OOB, kernel_stats.target_oob_cells);
		stats.inc(Statistics::ANCHORED_LANE_CELLS_QUERY_OOB, kernel_stats.query_oob_cells);
		stats.inc(Statistics::ANCHORED_LANE_CELLS_OUT_OF_BAND, kernel_stats.out_of_band_cells);
#endif
		stats.inc(Statistics::SWIPE_TASKS_TOTAL);
		stats.inc(Statistics::TIME_SW, timer.microseconds());
	}
	stats.inc(Statistics::TIME_ANCHORED_SWIPE, total.microseconds());
}

}

DISPATCH_6(void, anchored_swipe, const Extension* const*, begin, const Extension* const*, end, Loc, band_max, const int8_t*, common_scores, Statistics&, stats, const ExtensionCallback&, callback)

}

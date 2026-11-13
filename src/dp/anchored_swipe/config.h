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
#include <algorithm>
#include "basic/sequence.h"
#include "util/geo/geo.h"
#include "../long_profile/score_profile.h"
#include "align/new/extension.h"

namespace DP { namespace AnchoredSwipe {

struct Config {
	Score score_hint;
	Statistics& stats;
	Extension::Mode extension_mode;
	// The extensions are computed in the diagonal ranges of their anchors as set by the caller, who widens them by the
	// desired band. band_max has to be at least the band of every target, i.e. of the diagonal range of an extension
	// clipped to the sub-matrix it extends into. The anchored swipe sizes its data structures and profile paddings for
	// it. The caller sets it as a conservative bound over any extension that may be submitted to the queue (not only the
	// ones passed), so that extensions not known upon invocation can be computed within the same limits.
	Loc band_max;
};

/* If true, the kernel does not read the match scores from query profiles. Instead, it gets the scores of a
   target letter along the query by loading the letter's row of the score matrix and looking up the query
   letters in it with shuffle instructions, so that no profiles have to be built. */
static constexpr bool MATRIX_ROW_SCORES = true;

/* In MATRIX_ROW_SCORES mode, the query letters are framed by this letter, which scores PADDING_SCORE against
   every target letter in the score table (like the padding of the profiles). */
static constexpr Letter PADDING_LETTER = 30;
static constexpr int8_t PADDING_SCORE = -1;

/* If true, the kernel masks the cells of a channel that lie outside the band of its target, i.e. the rows of the kernel
   band at or beyond Target::band(), so that they do not contribute to the target's score (neither directly nor through
   gaps opened in them). The kernel band is the widest band of the targets loaded so far, rounded up to the channel
   count. Without the mask, a target is computed over that wider band, so that its result depends on which targets
   share the kernel with it. */
static constexpr bool MASK_TARGET_BAND = true;

/* Writes the score table for MATRIX_ROW_SCORES mode (see Options::score_table) of a score matrix given row-wise by
   target letter with 32 columns, like ScoreMatrix::matrix8() or Stats::TargetMatrix::scores, to out (32x32).
   Only the first `rows` rows are read from the matrix, the remaining rows of the table are zero. */
inline void make_score_table(const int8_t* matrix, int rows, int8_t* out) {
	std::copy(matrix, matrix + rows * 32, out);
	std::fill(out + rows * 32, out + 32 * 32, (int8_t)0);
	for (int l = 0; l < 32; ++l)
		out[(l << 5) + PADDING_LETTER] = PADDING_SCORE;
}

struct Options {
	const int16_t* const* profile, * const* profile_rev;
	// Score matrix for MATRIX_ROW_SCORES mode: 32x32, row-wise by target letter, with the PADDING_LETTER column set
	// (see make_score_table). Used for the targets that do not carry a score table of their own.
	const int8_t* score_table;
};

struct Stats {
	Stats():
		gross_cells(0),
		net_cells(0)
	{}
	int64_t gross_cells, net_cells;
#ifdef ANCHORED_SWIPE_STAT
	/* Lane cells computed by the kernel (kernel band x channels per column), each classified by the first of these
	   that applies: the channel is idle, the target coordinate is out of bounds, the query coordinate is out of
	   bounds, the cell lies outside the band of the channel's target. The remaining cells are useful. */
	int64_t lane_cells = 0, idle_cells = 0, target_oob_cells = 0, query_oob_cells = 0, out_of_band_cells = 0;
#endif
};

template<typename Score>
struct Target {
	Target()
	{}
	//Target(Sequence seq, Loc d_begin, Loc d_end, const Score* const* profile, Loc query_len, int64_t target_idx, bool reverse):
	Target(Sequence seq, Loc d_begin, Loc d_end, Loc query_start, Loc query_len, int64_t target_idx, bool reverse) :
		seq(seq),
		d_begin(d_begin),
		d_end(d_end),
		query_start(query_start),
		query_length(query_len),
		target_idx(target_idx),
		reverse(reverse),
		score(0),
		query_end(0),
		target_end(0)
	{}
	Sequence seq;
	Loc d_begin, d_end;
	//std::array<const Score* const*, 2> profile;
	Loc query_start, query_length;
	int64_t target_idx;
	bool reverse;
	Score score;
	Loc query_end, target_end;
	const LongScoreProfile<int16_t>* profile, *profile_rev;
	// MATRIX_ROW_SCORES mode: first letter of the padded query and of the padded reversed query.
	const Letter* query, *query_rev;
	// MATRIX_ROW_SCORES mode: score table the target is computed with, e.g. a compositionally adjusted matrix, in the
	// layout of Options::score_table. nullptr selects Options::score_table.
	const int8_t* score_table = nullptr;
	bool blank() const {
		return seq.length() == 0;
	}
	void reset() {
		seq = Sequence();
	}
	Loc band() const {
		return d_end - d_begin;
	}
	bool operator<(const Target& t) const {
		return band() < t.band();
	}
	static bool cmp_target_idx(const Target& a, const Target& b) {
		return a.target_idx < b.target_idx;
	}
	std::pair<int64_t, int64_t> cells() const {
		int64_t n = 0, g = 0;
		for (int j = 0; j < seq.length(); ++j) {
			int i0 = std::max(Geo::i(j, d_begin), 0),
				i1 = std::min(Geo::i(j, d_end), query_length);
			//assert(i1 - i0 >= 0);
			n += std::max(i1 - i0, 0);
			g += d_end - d_begin;
		}
		return { g,n };
	}
	int64_t gross_cells() const {
		return int64_t(d_end - d_begin) * seq.length();
	}
};

void anchored_swipe(std::vector<ExtensionPipeline::Extension>& batch, ExtensionPipeline::ExtensionQueue& queue, const DP::AnchoredSwipe::Config& cfg,
	const ExtensionPipeline::ExtensionCallback& callback);

}}
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
#include "basic/sequence.h"
#include "util/geo/geo.h"
#include "../long_profile/score_profile.h"

namespace DP { namespace AnchoredSwipe {

/* If true, the kernel does not read the match scores from query profiles. Instead, it gets the scores of a
   target letter along the query by loading the letter's row of the score matrix and looking up the query
   letters in it with shuffle instructions, so that no profiles have to be built. */
static constexpr bool MATRIX_ROW_SCORES = true;

/* In MATRIX_ROW_SCORES mode, the query letters are framed by this letter, which scores PADDING_SCORE against
   every target letter in the score table (like the padding of the profiles). */
static constexpr Letter PADDING_LETTER = 30;
static constexpr int8_t PADDING_SCORE = -1;

struct Options {
	const int16_t* const* profile, * const* profile_rev;
	// Score matrix for MATRIX_ROW_SCORES mode: 32x32, row-wise by target letter, with the PADDING_LETTER column set.
	const int8_t* score_table;
};

struct Stats {
	Stats():
		gross_cells(0),
		net_cells(0)
	{}
	int64_t gross_cells, net_cells;
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

}}
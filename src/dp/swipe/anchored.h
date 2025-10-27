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
#include <array>
#include "basic/sequence.h"
#include "../score_vector.h"
#include "util/simd/transpose.h"
#include "banded_matrix.h"
#include "config.h"
#include "util/geo/geo.h"
#include "util/util.h"
#include "util/data_structures/array.h"
#include "../score_vector_int8.h"
#include "../score_vector_int16.h"

using std::copy;
using std::array;
using std::numeric_limits;
using std::vector;
using std::pair;
using std::tie;

namespace DP { namespace AnchoredSwipe {

#if ARCH_ID == 2
	
namespace DISPATCH_ARCH {

static char blank[64] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
static constexpr Loc L = 13;

template<typename ScoreVector>
struct TargetIterator {
	enum { CHANNELS = ::DISPATCH_ARCH::ScoreTraits<ScoreVector>::CHANNELS };
	using Score = typename ::DISPATCH_ARCH::ScoreTraits<ScoreVector>::Score;
	TargetIterator(Target<Score>* targets, int64_t target_count, Loc target_len_max, DP::BandedSwipe::DISPATCH_ARCH::Matrix<ScoreVector>& matrix, const Options& options) :
		options(options),
		begin(targets),
		next(targets),
		end(targets + target_count),
		active(0),
		band(0)
	{
		while (active < CHANNELS && next < end) {
			int i = active;
			target_seqs[i] = Array<Letter>(target_len_max + 32 + 1);
			init_target(i);
			matrix.init_channel_diag(i, -Geo::i(0, targets[i].d_begin));
		}
		for (int i = active; i < CHANNELS; ++i)
			reset_channel(i);
		band = round_up(band, (Loc)CHANNELS);
		//matrix.init_channels_nw(-Geo::i(0, targets[0].d_begin), score_matrix.gap_open(), score_matrix.gap_extend());
	}
	inline bool init_target(int channel) {
		while (next->band() <= 0 && next < end) ++next;
		if (next == end)
			return false;
		target_idx[channel] = int(next - begin);
		targets[channel] = *next++;
		loc[channel] = 0;
		target_seqs[channel].assign(MASK_LETTER);
		if (targets[channel].reverse)
			target_seqs[channel].push_back_reversed(targets[channel].seq.data(), targets[channel].seq.end());
		else
			target_seqs[channel].push_back(targets[channel].seq.data(), targets[channel].seq.end());
		target_seqs[channel].push_back(32, MASK_LETTER);
		if (options.profile == nullptr) {
			if (targets[channel].reverse)
				for (int j = 0; j < AMINO_ACID_COUNT; ++j)
					profile_ptrs[channel][j] = targets[channel].profile_rev->get((int)j, targets[channel].query_start + Geo::i(0, targets[channel].d_begin) - 1);
			else
				for (int j = 0; j < AMINO_ACID_COUNT; ++j)
					profile_ptrs[channel][j] = targets[channel].profile->get((int)j, targets[channel].query_start + Geo::i(0, targets[channel].d_begin) - 1);
		}
		else {
			for (int j = 0; j < AMINO_ACID_COUNT; ++j)
				profile_ptrs[channel][j] = (const Score*)(targets[channel].reverse ? options.profile_rev[j] : options.profile[j])
				+ targets[channel].query_start + Geo::i(0, targets[channel].d_begin) - 1;
		}
		++active;
		band = std::max(band, targets[channel].band());
		return true;
	}
	inline void init_target_matrix(int channel, DP::BandedSwipe::DISPATCH_ARCH::Matrix<ScoreVector>& matrix, ScoreVector& max_score, ScoreVector& col_counter, ScoreVector& max_j) {
		matrix.init_channel_nw(channel, -Geo::i(0, targets[channel].d_begin), score_matrix.gap_open(), score_matrix.gap_extend());
		set_channel(max_score, channel, -1);
		set_channel(col_counter, channel, 0);
		set_channel(max_j, channel, -1);
	}
	inline void reset_channel(int channel) {
		if (options.profile == nullptr) {
			for (int j = 0; j < AMINO_ACID_COUNT; ++j)
				profile_ptrs[channel][j] = (const Score*)blank;
		}
		else
			for (int j = 0; j < AMINO_ACID_COUNT; ++j)
				profile_ptrs[channel][j] = (const Score*)options.profile[0];
	}
	inline void next_block(DP::BandedSwipe::DISPATCH_ARCH::Matrix<ScoreVector>& matrix, ScoreVector& max_score, ScoreVector& max_i, ScoreVector& max_j, ScoreVector& col_counter) {
		for (int i = 0; i < CHANNELS; ++i) {
			if (targets[i].blank()) {
				std::fill(letters[i].begin(), letters[i].end(), MASK_LETTER);
				continue;
			}
			if (loc[i] >= targets[i].seq.length() + 1) {
				const ::Score score = max_score[i];
				if (score >= 0) {
					begin[target_idx[i]].score = score + 1;
					const Score j1 = max_j[i];
					if (j1 < numeric_limits<Score>::max()) {
						begin[target_idx[i]].target_end = (Loc)j1 + 1 - 1;
						begin[target_idx[i]].query_end = Geo::i((Loc)j1, targets[i].d_begin) + (Loc)max_i[i] + 1 - 1;
						assert(begin[target_idx[i]].target_end > 0 && begin[target_idx[i]].query_end > 0);
					}
					else
						begin[target_idx[i]].score = numeric_limits<Score>::max();
				}
				--active;
				targets[i].reset();
				if (next < end) {
					if (!init_target(i)) {
						std::fill(letters[i].begin(), letters[i].end(), MASK_LETTER);
						reset_channel(i);
						continue;
					}
					init_target_matrix(i, matrix, max_score, col_counter, max_j);
					band = round_up(band, (Loc)CHANNELS);
				}
				else {
					std::fill(letters[i].begin(), letters[i].end(), MASK_LETTER);
					reset_channel(i);
					continue;
				}
			}
			copy(target_seqs[i].data() + loc[i], target_seqs[i].data() + loc[i] + L, letters[i].data());
			loc[i] += L;
			if (profile_ptrs[i][0] != (const Score*)blank)
				for (int j = 0; j < AMINO_ACID_COUNT; ++j)
					profile_ptrs[i][j] += L;
		}
	}
	inline array<const Score*, CHANNELS> column_ptrs(int k) {
		array<const Score*, CHANNELS> prof_ptr;
		for (int i = 0; i < CHANNELS; ++i) {
			if(profile_ptrs[i][0] == (const Score*)blank) {
				prof_ptr[i] = (const Score*)blank;
				continue;
			}
			const Letter l = letter_mask(letters[i][k + L]);
			prof_ptr[i] = profile_ptrs[i][(int)l] + k;
		}
		return prof_ptr;
	}
	size_t net_cells(int k) const {
		size_t n = 0;
		for (int i = 0; i < CHANNELS; ++i)
			if (!targets[i].blank() && (loc[i] + k < targets[i].seq.length())) {
				int j = loc[i] + k, i0 = std::max(Geo::i(j, targets[i].d_begin), 0),
					i1 = std::min(Geo::i(j, targets[i].d_end), targets[i].query_length);
				//assert(i1 - i0 >= 0);
				n += std::max(i1 - i0, 0);
			}
		return n;
	}
	const Options& options;
	array<Target<Score>, CHANNELS> targets;
	array<Array<Letter>, CHANNELS> target_seqs;
	Target<Score>* begin, * next, * end;
	int active;
	array<array<const Score*, AMINO_ACID_COUNT>, CHANNELS> profile_ptrs;
	array<Loc, CHANNELS> loc;
	array<array<Letter, L>, CHANNELS> letters;
	array<char, 8192> padding;
	array<int, CHANNELS> target_idx;
	Loc band;
};

template<typename ScoreVector>
Stats FLATTEN smith_waterman(DP::AnchoredSwipe::Target<typename ::DISPATCH_ARCH::ScoreTraits<ScoreVector>::Score>* targets, int64_t target_count, const Options& options) {
	using Score = typename ::DISPATCH_ARCH::ScoreTraits<ScoreVector>::Score;
	const Loc CHANNELS = ::DISPATCH_ARCH::ScoreTraits<ScoreVector>::CHANNELS;
	constexpr Score SCORE_MIN = numeric_limits<Score>::min();
	if (target_count == 0)
		return Stats();

	alignas(32) Score scores[CHANNELS * CHANNELS];
	Loc band_max, target_len_max;
	tie(band_max, target_len_max) = limits(targets, target_count);
	DP::BandedSwipe::DISPATCH_ARCH::Matrix<ScoreVector> matrix(round_up(band_max, CHANNELS), 0, ScoreVector(SCORE_MIN));
	assert(round_up(band_max, CHANNELS) <= numeric_limits<Score>::max());
	TargetIterator<ScoreVector> target_it(targets, target_count, target_len_max, matrix, options);
	const ScoreVector go = ScoreVector(score_matrix.gap_open() + score_matrix.gap_extend()),
		ge = ScoreVector(score_matrix.gap_extend()), one = ScoreVector(1);
	ScoreVector max_score(-1), col_counter(0), max_j(-1), max_i(0);
	Stats stats;

	while(target_it.next_block(matrix, max_score, max_i, max_j, col_counter), target_it.active > 0) {
		const int band = target_it.band;
		for (int k = -L; k < 0; ++k) {
#ifdef DP_STAT
			stats.gross_cells += (size_t)band * CHANNELS;
			stats.net_cells += target_it.net_cells(k);
#endif

			typename DP::BandedSwipe::DISPATCH_ARCH::Matrix<ScoreVector>::ColumnIterator it(matrix.begin(0, 0));
			array<const Score*, CHANNELS> prof_ptr = target_it.column_ptrs(k);
			ScoreVector vgap = ScoreVector(SCORE_MIN), hgap = ScoreVector(), col_best = ScoreVector(SCORE_MIN), row_counter(0), col_max_i(0);

			for (int i = 0; i < band;) {
				transpose_offset(prof_ptr.data(), CHANNELS, i / CHANNELS, scores, __m256i());
				const Score* score_ptr = scores;

				do {
					hgap = it.hgap();
					ScoreVector match_scores(score_ptr);
					ScoreVector score = it.diag() + match_scores;
					score = max(score, hgap);
					score = max(score, vgap);
					ScoreVector open = score - go;
					const ScoreVector gt_mask = score > col_best;
					col_max_i = blend(col_max_i, row_counter, gt_mask);
					row_counter += one;
					col_best = max(col_best, score);
					vgap -= ge;
					hgap -= ge;
					vgap = max(vgap, open);
					hgap = max(hgap, open);
					it.set_hgap(hgap);
					it.set_score(score);
					++it;
					score_ptr += CHANNELS;
					++i;
				} while ((i & (CHANNELS - 1)) != 0);
			}
			const ScoreVector gt_mask = col_best > max_score;
			max_j = blend(max_j, col_counter, gt_mask);
			max_i = blend(max_i, col_max_i, gt_mask);
			max_score = max(max_score, col_best);
			col_counter += one;
		}
	}
	return stats;
}

}

#endif

}}
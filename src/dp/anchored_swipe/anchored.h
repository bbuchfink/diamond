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
#include <array>
#include <memory>
#include <type_traits>
#include "basic/sequence.h"
#include "../score_vector.h"
#include "util/simd/transpose.h"
#include "../swipe/banded_matrix.h"
#include "config.h"
#include "util/geo/geo.h"
#include "util/util.h"
#include "../score_vector_int8.h"
#include "../score_vector_int16.h"

using std::copy;
using std::array;
using std::numeric_limits;
using std::vector;
using std::pair;
using std::tie;

namespace DP { namespace AnchoredSwipe {

/* The kernel runs on the 16 bit score vector of whatever architecture variant this
   translation unit is compiled into: 16 channels on AVX2/AVX512, 8 on SSE2/SSE4.1
   and NEON, and the scalar single lane fallback everywhere else. */
	
namespace DISPATCH_ARCH {

static constexpr Loc L = 13;

/* MATRIX_ROW_SCORES mode: writes the scores of the letters seq[0..CHANNELS) against the score matrix row `row`
   (32 entries indexed by letter) to out[0..CHANNELS), sign extended to 16 bit. Reads exactly CHANNELS letters.
   The overloads are selected by the register type of the 16 bit score vector, like transpose_offset. */

static inline void score_row(const int8_t* row, const Letter* seq, int16_t* out, int16_t) {
	*out = row[(int)letter_mask(*seq)];
}

#ifdef __AVX2__

static inline void score_row(const int8_t* row, const Letter* seq, int16_t* out, __m256i) {
	const __m128i s = _mm_and_si128(_mm_loadu_si128((const __m128i*)seq), _mm_set1_epi8(LETTER_MASK));
	const __m128i high_mask = _mm_slli_epi16(_mm_and_si128(s, _mm_set1_epi8('\x10')), 3);
	const __m128i s1 = _mm_shuffle_epi8(_mm_loadu_si128((const __m128i*)row), _mm_or_si128(s, high_mask));
	const __m128i s2 = _mm_shuffle_epi8(_mm_loadu_si128((const __m128i*)(row + 16)), _mm_or_si128(s, _mm_xor_si128(high_mask, _mm_set1_epi8('\x80'))));
	_mm256_storeu_si256((__m256i*)out, _mm256_cvtepi8_epi16(_mm_or_si128(s1, s2)));
}

#endif

#ifdef __SSE2__

static inline void score_row(const int8_t* row, const Letter* seq, int16_t* out, __m128i) {
#ifdef __SSSE3__
	const __m128i s = _mm_and_si128(_mm_loadl_epi64((const __m128i*)seq), _mm_set1_epi8(LETTER_MASK));
	const __m128i high_mask = _mm_slli_epi16(_mm_and_si128(s, _mm_set1_epi8('\x10')), 3);
	const __m128i s1 = _mm_shuffle_epi8(_mm_loadu_si128((const __m128i*)row), _mm_or_si128(s, high_mask));
	const __m128i s2 = _mm_shuffle_epi8(_mm_loadu_si128((const __m128i*)(row + 16)), _mm_or_si128(s, _mm_xor_si128(high_mask, _mm_set1_epi8('\x80'))));
	const __m128i r = _mm_or_si128(s1, s2);
#ifdef __SSE4_1__
	_mm_storeu_si128((__m128i*)out, _mm_cvtepi8_epi16(r));
#else
	_mm_storeu_si128((__m128i*)out, _mm_unpacklo_epi8(r, _mm_cmpgt_epi8(_mm_setzero_si128(), r)));
#endif
#else
	for (int i = 0; i < 8; ++i)
		out[i] = row[(int)letter_mask(seq[i])];
#endif
}

#endif

#ifdef __ARM_NEON

static inline void score_row(const int8_t* row, const Letter* seq, int16_t* out, int16x8_t) {
	const int8x8_t s = vand_s8(vld1_s8(seq), vdup_n_s8(LETTER_MASK));
#ifdef __aarch64__
	int8x16x2_t t;
	t.val[0] = vld1q_s8(row);
	t.val[1] = vld1q_s8(row + 16);
	const int8x8_t r = vqtbl2_s8(t, vreinterpret_u8_s8(s));
#else
	int8x8x4_t t;
	for (int i = 0; i < 4; ++i)
		t.val[i] = vld1_s8(row + 8 * i);
	const int8x8_t r = vtbl4_s8(t, s);
#endif
	vst1q_s16(out, vmovl_s8(r));
}

#endif

// Target sequence of a channel, grown as needed. Same size as the fixed Array it replaces, which keeps the layout of TargetIterator.
struct TargetSeqBuffer {
	std::unique_ptr<Letter[]> data;
	Loc capacity = 0;
};

template<typename ScoreVector>
struct TargetIterator {
	enum { CHANNELS = ::DISPATCH_ARCH::ScoreTraits<ScoreVector>::CHANNELS };
	using Score = typename ::DISPATCH_ARCH::ScoreTraits<ScoreVector>::Score;
	TargetIterator(Target<Score>* targets, int64_t target_count, DP::BandedSwipe::DISPATCH_ARCH::Matrix<ScoreVector>& matrix, const Options& options) :
		options(options),
		begin(targets),
		next(targets),
		end(targets + target_count),
		active(0),
		band(0),
		blank_profile(std::max(matrix.band(), (int)CHANNELS), (Score)0),
		blank_query(std::max(matrix.band(), (int)CHANNELS), (Letter)0)
	{
		assert(!MATRIX_ROW_SCORES || (options.profile == nullptr && options.score_table != nullptr));
		zero_row.fill(0);
		while (active < CHANNELS && next < end) {
			int i = active;
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
		TargetSeqBuffer& buf = target_seqs[channel];
		const Sequence& target_seq = targets[channel].seq;
		const Loc size = target_seq.length() + 1 + L;
		if (size > buf.capacity) {
			buf.capacity = size * 2;
			buf.data.reset(new Letter[buf.capacity]);
		}
		Letter* seq = buf.data.get();
		*seq++ = MASK_LETTER;
		if (targets[channel].reverse)
			std::reverse_copy(target_seq.data(), target_seq.end(), seq);
		else
			std::copy(target_seq.data(), target_seq.end(), seq);
		std::fill(seq + target_seq.length(), seq + target_seq.length() + L, MASK_LETTER);
		if (MATRIX_ROW_SCORES)
			query_ptrs[channel] = (targets[channel].reverse ? targets[channel].query_rev : targets[channel].query)
				+ targets[channel].query_start + Geo::i(0, targets[channel].d_begin) - 1;
		else if (options.profile == nullptr) {
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
		if (MATRIX_ROW_SCORES)
			query_ptrs[channel] = blank_query.data();
		else if (options.profile == nullptr) {
			for (int j = 0; j < AMINO_ACID_COUNT; ++j)
				profile_ptrs[channel][j] = blank_profile.data();
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
			copy(target_seqs[i].data.get() + loc[i], target_seqs[i].data.get() + loc[i] + L, letters[i].data());
			loc[i] += L;
			if (MATRIX_ROW_SCORES) {
				if (query_ptrs[i] != blank_query.data())
					query_ptrs[i] += L;
			}
			else if (profile_ptrs[i][0] != blank_profile.data())
				for (int j = 0; j < AMINO_ACID_COUNT; ++j)
					profile_ptrs[i][j] += L;
		}
	}
	inline array<const Score*, CHANNELS> column_ptrs(int k) {
		array<const Score*, CHANNELS> prof_ptr;
		for (int i = 0; i < CHANNELS; ++i) {
			if(profile_ptrs[i][0] == blank_profile.data()) {
				prof_ptr[i] = blank_profile.data();
				continue;
			}
			const Letter l = letter_mask(letters[i][k + L]);
			prof_ptr[i] = profile_ptrs[i][(int)l] + k;
		}
		return prof_ptr;
	}
	// MATRIX_ROW_SCORES mode: the score matrix rows of the target letters of column k, and the query letters of row 0.
	inline void column_rows(int k, array<const int8_t*, CHANNELS>& rows, array<const Letter*, CHANNELS>& query) {
		for (int i = 0; i < CHANNELS; ++i) {
			if (query_ptrs[i] == blank_query.data()) {
				rows[i] = zero_row.data();
				query[i] = blank_query.data();
				continue;
			}
			const Letter l = letter_mask(letters[i][k + L]);
			rows[i] = options.score_table + ((int)l << 5);
			query[i] = query_ptrs[i] + k;
		}
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
	array<TargetSeqBuffer, CHANNELS> target_seqs;
	Target<Score>* begin, * next, * end;
	int active;
	array<array<const Score*, AMINO_ACID_COUNT>, CHANNELS> profile_ptrs;
	array<const Letter*, CHANNELS> query_ptrs;
	array<Loc, CHANNELS> loc;
	array<array<Letter, L>, CHANNELS> letters;
	array<char, 8192> padding;
	array<int, CHANNELS> target_idx;
	Loc band;
	vector<Score> blank_profile;
	vector<Letter> blank_query;
	array<int8_t, 32> zero_row;
};

/* Profile padding for which the kernel does not read outside the profiles of targets with band() <= band_max
   and d_end >= 1, i.e. whose band contains the diagonal of the start cell. In MATRIX_ROW_SCORES mode, the kernel
   reads the query letters at the same positions, so the padded queries need the same padding.
   For a target, the kernel reads the profile at query position query_start + d_begin - 1 + c + i for the
   columns c in [0, round_up(tlen + 1, L)) and the rows i in [0, band), where band is the widest band of the
   targets loaded so far rounded up to the channel count, since every channel is computed over the full kernel
   band. As the target is clipped to tlen <= query_length - d_begin, the read beyond the query end is at most
   round_up(band_max, CHANNELS) + L - 2 positions. The read before the query start is at most -d_begin + 1 <= band_max positions. */
template<typename ScoreVector>
int64_t profile_padding(Loc band_max) {
	return round_up(band_max, (Loc)::DISPATCH_ARCH::ScoreTraits<ScoreVector>::CHANNELS) + L - 2;
}

/* Computes the targets, which have to satisfy band() <= band_max and d_end >= 1. The profiles have to be padded by
   at least profile_padding(band_max). */
template<typename ScoreVector>
Stats FLATTEN smith_waterman(DP::AnchoredSwipe::Target<typename ::DISPATCH_ARCH::ScoreTraits<ScoreVector>::Score>* targets, int64_t target_count, Loc band_max, const Options& options) {
	using Score = typename ::DISPATCH_ARCH::ScoreTraits<ScoreVector>::Score;
	const Loc CHANNELS = ::DISPATCH_ARCH::ScoreTraits<ScoreVector>::CHANNELS;
	constexpr Score SCORE_MIN = numeric_limits<Score>::min();
	if (target_count == 0)
		return Stats();

	alignas(32) Score scores[CHANNELS * CHANNELS];
	// MATRIX_ROW_SCORES mode: scores of each channel along the query, transposed into scores.
	alignas(32) Score row_scores[CHANNELS * CHANNELS];
	array<const Score*, CHANNELS> row_score_ptrs;
	for (int c = 0; c < CHANNELS; ++c)
		row_score_ptrs[c] = row_scores + c * CHANNELS;
	static_assert(!MATRIX_ROW_SCORES || std::is_same<Score, int16_t>::value, "MATRIX_ROW_SCORES requires 16 bit scores");
	DP::BandedSwipe::DISPATCH_ARCH::Matrix<ScoreVector> matrix(round_up(band_max, CHANNELS), 0, nullptr, ScoreVector(SCORE_MIN));
	assert(round_up(band_max, CHANNELS) <= numeric_limits<Score>::max());
	TargetIterator<ScoreVector> target_it(targets, target_count, matrix, options);
	const ScoreVector go = ScoreVector(score_matrix.gap_open() + score_matrix.gap_extend()),
		ge = ScoreVector(score_matrix.gap_extend()), one = ScoreVector(1);
	ScoreVector max_score(-1), col_counter(0), max_j(-1), max_i(0);
	Stats stats;

	while(target_it.next_block(matrix, max_score, max_i, max_j, col_counter), target_it.active > 0) {
		const int band = target_it.band;
		assert(band <= (int)target_it.blank_profile.size());
		for (int k = -L; k < 0; ++k) {
#ifdef DP_STAT
			stats.gross_cells += (size_t)band * CHANNELS;
			stats.net_cells += target_it.net_cells(k);
#endif

			typename DP::BandedSwipe::DISPATCH_ARCH::Matrix<ScoreVector>::ColumnIterator it(matrix.begin(0, 0));
			array<const Score*, CHANNELS> prof_ptr;
			array<const int8_t*, CHANNELS> rows;
			array<const Letter*, CHANNELS> query;
			if (MATRIX_ROW_SCORES)
				target_it.column_rows(k, rows, query);
			else
				prof_ptr = target_it.column_ptrs(k);
			ScoreVector vgap = ScoreVector(SCORE_MIN), hgap = ScoreVector(), col_best = ScoreVector(SCORE_MIN), row_counter(0), col_max_i(0);

			for (int i = 0; i < band;) {
				if (MATRIX_ROW_SCORES) {
					for (int c = 0; c < CHANNELS; ++c)
						score_row(rows[c], query[c] + i, row_scores + c * CHANNELS, typename ScoreVector::Register());
					transpose_offset(row_score_ptrs.data(), CHANNELS, 0, scores, typename ScoreVector::Register());
				}
				else
					transpose_offset(prof_ptr.data(), CHANNELS, i / CHANNELS, scores, typename ScoreVector::Register());
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

}}
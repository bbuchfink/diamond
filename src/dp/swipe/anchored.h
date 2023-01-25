#pragma once
#include <array>
#include "../../basic/sequence.h"
#include "../score_vector.h"
#include "../../util/simd/transpose32x32.h"
#include "../../util/simd/transpose16x16.h"
#include "banded_matrix.h"
#include "cell_update.h"
#include "config.h"
#include "../../util/geo/geo.h"
#include "../score_vector_int8.h"
#include "../score_vector_int16.h"
#include "../../util/util.h"
#include "../../util/data_structures/array.h"

using std::copy;
using std::array;
using std::numeric_limits;
using std::vector;
using std::pair;
using std::tie;

namespace DP { namespace AnchoredSwipe {

struct Options {
	const int16_t* const* profile, * const* profile_rev;
};

#if ARCH_ID == 2
	
namespace DISPATCH_ARCH {

static char blank[64] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
static constexpr Loc L = 13;

template<typename Score>
pair<Loc, Loc> limits(const Target<Score>* targets, size_t count) {
	int band = 0, target_len = 0;
	for (const Target<Score>* i = targets; i < targets + count; ++i) {
		assert(i->band() > 0);
		band = std::max(band, i->band());
		target_len = std::max(target_len, i->seq.length());
	}
	return { band,target_len };
}

template<typename Sv>
struct TargetIterator {
	enum { CHANNELS = ::DISPATCH_ARCH::ScoreTraits<Sv>::CHANNELS };
	using Score = typename ::DISPATCH_ARCH::ScoreTraits<Sv>::Score;
	TargetIterator(Target<Score>* targets, int64_t target_count, Loc target_len_max, DP::BandedSwipe::DISPATCH_ARCH::Matrix<Sv>& matrix, const Options& options) :
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
		for (int j = 0; j < AMINO_ACID_COUNT; ++j)
			profile_ptrs[channel][j] = (const Score*)(targets[channel].reverse ? options.profile_rev[j] : options.profile[j])
				+ targets[channel].query_start + Geo::i(0, targets[channel].d_begin) - 1;		
		++active;
		band = std::max(band, targets[channel].band());
		return true;
	}
	inline void init_target_matrix(int channel, DP::BandedSwipe::DISPATCH_ARCH::Matrix<Sv>& matrix, Sv& max_score, Sv& col_counter, Sv& max_j) {
		matrix.init_channel_nw(channel, -Geo::i(0, targets[channel].d_begin), score_matrix.gap_open(), score_matrix.gap_extend());
		set_channel(max_score, channel, -1);
		set_channel(col_counter, channel, 0);
		set_channel(max_j, channel, -1);
	}
	inline void reset_channel(int channel) {
		for (int j = 0; j < AMINO_ACID_COUNT; ++j)
			profile_ptrs[channel][j] = (const Score*)options.profile[0];
	}
	inline void next_block(DP::BandedSwipe::DISPATCH_ARCH::Matrix<Sv>& matrix, Sv& max_score, Sv& max_i, Sv& max_j, Sv& col_counter) {
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
			for (int j = 0; j < AMINO_ACID_COUNT; ++j)
				profile_ptrs[i][j] += L;
		}
	}
	inline array<const Score*, CHANNELS> column_ptrs(int k) {
		array<const Score*, CHANNELS> prof_ptr;
		for (int i = 0; i < CHANNELS; ++i) {
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

template<typename Sv>
Stats FLATTEN smith_waterman(DP::AnchoredSwipe::Target<typename ::DISPATCH_ARCH::ScoreTraits<Sv>::Score>* targets, int64_t target_count, const Options& options) {
	using Score = typename ::DISPATCH_ARCH::ScoreTraits<Sv>::Score;
	const Loc CHANNELS = ::DISPATCH_ARCH::ScoreTraits<Sv>::CHANNELS;
	constexpr Score SCORE_MIN = numeric_limits<Score>::min();
	if (target_count == 0)
		return Stats();

	alignas(32) Score scores[CHANNELS * CHANNELS];
	Loc band_max, target_len_max;
	tie(band_max, target_len_max) = limits(targets, target_count);
	DP::BandedSwipe::DISPATCH_ARCH::Matrix<Sv> matrix(round_up(band_max, CHANNELS), 0, Sv(SCORE_MIN));
	assert(round_up(band_max, CHANNELS) <= numeric_limits<Score>::max());
	TargetIterator<Sv> target_it(targets, target_count, target_len_max, matrix, options);
	const Sv go = Sv(score_matrix.gap_open() + score_matrix.gap_extend()),
		ge = Sv(score_matrix.gap_extend()), one = Sv(1);
	Sv max_score(-1), col_counter(0), max_j(-1), max_i(0);
	Stats stats;

	while(target_it.next_block(matrix, max_score, max_i, max_j, col_counter), target_it.active > 0) {
		const int band = target_it.band;
		for (int k = -L; k < 0; ++k) {
#ifdef DP_STAT
			stats.gross_cells += (size_t)band * CHANNELS;
			stats.net_cells += target_it.net_cells(k);
#endif

			typename DP::BandedSwipe::DISPATCH_ARCH::Matrix<Sv>::ColumnIterator it(matrix.begin(0, 0));
			array<const Score*, CHANNELS> prof_ptr = target_it.column_ptrs(k);
			Sv vgap = Sv(SCORE_MIN), hgap = Sv(), col_best = Sv(SCORE_MIN), row_counter(0), col_max_i(0);

			for (int i = 0; i < band;) {
			    transpose_offset(prof_ptr.data(), CHANNELS, i/CHANNELS, scores, __m256i());
				const Score* score_ptr = scores;

				do {
					hgap = it.hgap();
					Sv match_scores(score_ptr);
					Sv score = it.diag() + match_scores;
					score = max(score, hgap);
					score = max(score, vgap);
					Sv open = score - go;
					const Sv gt_mask = score > col_best;
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
			const Sv gt_mask = col_best > max_score;
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


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
#include <array>
#include <limits>
#include "score_profile.h"
#include "../score_vector.h"
#include "../score_vector_int8.h"
#include "../score_vector_int16.h"
#include "util/simd/dispatch.h"
#include "stats/score_matrix.h"

using std::array;
using std::vector;

namespace DP { namespace DISPATCH_ARCH {

static constexpr int PADDING_SCORE = -1;

#if defined(__SSE4_1__) | defined(__ARM_NEON)

using Sv8 = ::DISPATCH_ARCH::ScoreVector<int8_t, 0>;
using Sv16 = ::DISPATCH_ARCH::ScoreVector<int16_t, 0>;

// Number of query positions scored per vector lookup.
static constexpr Loc CHUNK = ::DISPATCH_ARCH::ScoreTraits<Sv8>::CHANNELS;

// Writes the scores of the CHUNK letters at seq against one row of the score table, plus the composition bias if cbs is set.
static inline void score_chunk(const Sv8::Table& row, const Letter* seq, const int8_t* cbs, int8_t* out) {
	Sv8 s(row, seq);
	if (cbs)
		s += Sv8(cbs);
	s.store(out);
}

static inline void score_chunk(const Sv8::Table& row, const Letter* seq, const int8_t* cbs, int16_t* out) {
	store_expanded(Sv8(row, seq), out);
	if (!cbs)
		return;
	// The bias is added in 16 bit so that it can not saturate at the 8 bit range.
	constexpr Loc C16 = ::DISPATCH_ARCH::ScoreTraits<Sv16>::CHANNELS;
	for (Loc k = 0; k < CHUNK; k += C16)
		(Sv16(out + k) + Sv16::load_expanded(cbs + k)).store(out + k);
}

#endif

/* Builds the profile of seq against a 32x32 score table given row-wise by target letter.
   The composition bias cbs (optional) is added for the true amino acids only. */
template<typename Score>
static LongScoreProfile<Score> make_profile(Sequence seq, const int8_t* matrix, const int8_t* cbs, int64_t padding)
{
	LongScoreProfile<Score> p(padding);
	const Loc n = seq.length();
	const Letter* letters = seq.data();

#if defined(__SSE4_1__) | defined(__ARM_NEON)
	// The last partial chunk is scored from zero padded copies so that nothing is read beyond the sequence.
	constexpr Loc C16 = ::DISPATCH_ARCH::ScoreTraits<Sv16>::CHANNELS;
	static_assert(CHUNK % C16 == 0, "8 bit channel count must be a multiple of the 16 bit one");
	const Loc full = n - n % CHUNK;
	array<Letter, CHUNK> tail_seq;
	array<int8_t, CHUNK> tail_cbs;
	array<Score, CHUNK> tail_out;
	tail_seq.fill(0);
	tail_cbs.fill(0);
	std::copy(letters + full, letters + n, tail_seq.begin());
	if (cbs)
		std::copy(cbs + full, cbs + n, tail_cbs.begin());
#endif

	for (int l = 0; l < AMINO_ACID_COUNT; ++l) {
		vector<Score>& v = p.data[l];
		v.resize(n + 2 * p.padding);
		std::fill(v.begin(), v.begin() + p.padding, (Score)PADDING_SCORE);
		std::fill(v.end() - p.padding, v.end(), (Score)PADDING_SCORE);
		Score* out = v.data() + p.padding;
		const int8_t* row = &matrix[l << 5];
		const bool with_cbs = cbs && l < TRUE_AA;
#if defined(__SSE4_1__) | defined(__ARM_NEON)
		const Sv8::Table table(row);
		for (Loc i = 0; i < full; i += CHUNK)
			score_chunk(table, letters + i, with_cbs ? cbs + i : nullptr, out + i);
		if (full < n) {
			score_chunk(table, tail_seq.data(), with_cbs ? tail_cbs.data() : nullptr, tail_out.data());
			std::copy(tail_out.begin(), tail_out.begin() + (n - full), out + full);
		}
#else
		const int lo = std::numeric_limits<Score>::min(), hi = std::numeric_limits<Score>::max();
		for (Loc i = 0; i < n; ++i) {
			const int s = (int)row[(int)letter_mask(letters[i])] + (with_cbs ? (int)cbs[i] : 0);
			out[i] = (Score)std::min(std::max(s, lo), hi);
		}
#endif
	}
	return p;
}

LongScoreProfile<int8_t> make_profile8(Sequence seq, const int8_t* cbs, int64_t padding) {
	return make_profile<int8_t>(seq, score_matrix.matrix8(), cbs, padding);
}

LongScoreProfile<int16_t> make_profile16(Sequence seq, const int8_t* cbs, int64_t padding, const ::ScoreMatrix* matrix) {
	return make_profile<int16_t>(seq, matrix->matrix8(), cbs, padding);
}

LongScoreProfile<int16_t> make_profile16(Sequence seq, const Stats::TargetMatrix& matrix, int64_t padding) {
	return make_profile<int16_t>(seq, matrix.scores.data(), nullptr, padding);
}

}

DISPATCH_3(LongScoreProfile<int8_t>, make_profile8, Sequence, seq, const int8_t*, cbs, int64_t, padding);
DISPATCH_4(LongScoreProfile<int16_t>, make_profile16, Sequence, seq, const int8_t*, cbs, int64_t, padding, const ::ScoreMatrix*, matrix);
DISPATCH_3(LongScoreProfile<int16_t>, make_profile16, Sequence, seq, const Stats::TargetMatrix&, matrix, int64_t, padding);

}

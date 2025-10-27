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

#include "score_profile.h"
#include "score_vector.h"
#include "score_vector_int8.h"
#include "score_vector_int16.h"
#include "util/util.h"
#include "util/simd/dispatch.h"
#include "stats/score_matrix.h"

using std::array;

namespace DP { namespace DISPATCH_ARCH {

template<typename Score>
LongScoreProfile<Score> make_profile(Sequence seq, const int8_t* cbs, int64_t padding, const ScoreMatrix& matrix)
{
	LongScoreProfile<Score> p;
	p.padding = std::max(padding, (int64_t)LongScoreProfile<Score>::DEFAULT_PADDING);
	for (unsigned l = 0; l < AMINO_ACID_COUNT; ++l) {
		const int8_t* scores = &matrix.matrix8()[l << 5];
		p.data[l].reserve(round_up(seq.length(), 32) + 2 * p.padding);
		p.data[l].insert(p.data[l].end(), p.padding, -1);
#if ARCH_ID == 2
		using Sv = ::DISPATCH_ARCH::ScoreVector<int8_t, 0>;
		constexpr auto CHANNELS = ::DISPATCH_ARCH::ScoreTraits<Sv>::CHANNELS;
		alignas(32) array<Score, CHANNELS> buf;
		for (Loc i = 0; i < seq.length(); i += CHANNELS) {
			__m256i s = _mm256_loadu_si256((const __m256i*)(seq.data() + i));
			Sv scores(l, s, matrix.matrix8_low(), matrix.matrix8_high());
			if (cbs && l < TRUE_AA)
				scores += Sv(cbs + i);
			store_expanded(scores, buf.data());
			p.data[l].insert(p.data[l].end(), buf.begin(), buf.end());
		}
		p.data[l].erase(p.data[l].end() - round_up(seq.length(), (Loc)CHANNELS) + seq.length(), p.data[l].end());
#else
		for (Loc i = 0; i < seq.length(); ++i) {
			int8_t score = scores[(int)seq[i]];
			if (cbs)
				score += cbs[i];
			p.data[l].push_back(score);
		}
#endif
		p.data[l].insert(p.data[l].end(), p.padding, -1);
	}
	return p;
}

template<typename Score>
LongScoreProfile<Score> make_profile(Sequence seq, const Stats::TargetMatrix& matrix, int64_t padding)
{
	LongScoreProfile<Score> p;
	p.padding = std::max(padding, (int64_t)LongScoreProfile<Score>::DEFAULT_PADDING);
	for (int l = 0; l < AMINO_ACID_COUNT; ++l) {
		const int8_t* scores = &matrix.scores[l << 5];
		p.data[l].reserve(round_up(seq.length(), 32) + 2 * p.padding);
		p.data[l].insert(p.data[l].end(), p.padding, -1);
/*#if ARCH_ID == 2
		using Sv = ::DISPATCH_ARCH::ScoreVector<int8_t, 0>;
		constexpr auto CHANNELS = ::DISPATCH_ARCH::ScoreTraits<Sv>::CHANNELS;
		alignas(32) array<Score, CHANNELS> buf;
		const int8_t* low = matrix.scores_low.data, * high = matrix.scores_high.data;
		for (Loc i = 0; i < seq.length(); i += CHANNELS) {
			__m256i s = _mm256_loadu_si256((const __m256i*)(seq.data() + i));
			Sv scores(l, s, low, high);
			store_expanded(scores, buf.data());
			p.data[l].insert(p.data[l].end(), buf.begin(), buf.end());
		}
		p.data[l].erase(p.data[l].end() - round_up(seq.length(), (Loc)CHANNELS) + seq.length(), p.data[l].end());
#else*/
		for (Loc i = 0; i < seq.length(); ++i) {
			int8_t score = scores[(int)seq[i]];
			p.data[l].push_back(score);
		}
//#endif
		p.data[l].insert(p.data[l].end(), p.padding, -1);
	}
	return p;
}

LongScoreProfile<int8_t> make_profile8(Sequence seq, const int8_t* cbs, int64_t padding) {
	return make_profile<int8_t>(seq, cbs, padding, score_matrix);
}

LongScoreProfile<int16_t> make_profile16(Sequence seq, const int8_t* cbs, int64_t padding, const ::ScoreMatrix* matrix) {
	return make_profile<int16_t>(seq, cbs, padding, *matrix);
}

LongScoreProfile<int16_t> make_profile16(Sequence seq, const Stats::TargetMatrix& matrix, int64_t padding) {
	return make_profile<int16_t>(seq, matrix, padding);
}

}

DISPATCH_3(LongScoreProfile<int8_t>, make_profile8, Sequence, seq, const int8_t*, cbs, int64_t, padding);
DISPATCH_4(LongScoreProfile<int16_t>, make_profile16, Sequence, seq, const int8_t*, cbs, int64_t, padding, const ::ScoreMatrix*, matrix);
DISPATCH_3(LongScoreProfile<int16_t>, make_profile16, Sequence, seq, const Stats::TargetMatrix&, matrix, int64_t, padding);

}
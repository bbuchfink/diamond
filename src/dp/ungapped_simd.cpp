/****
DIAMOND protein aligner
Copyright (C) 2020 Max Planck Society for the Advancement of Science e.V.

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

#include <assert.h>
#include <algorithm>
#include "score_vector_int8.h"
#include "../util/simd/vector.h"
#include "ungapped_simd.h"
#include "../util/simd/transpose.h"
#include "ungapped.h"

using namespace DISPATCH_ARCH;

namespace DP { namespace DISPATCH_ARCH {

void window_ungapped(const Letter *query, const Letter **subjects, int subject_count, int window, int *out) {
#ifdef __SSE4_1__
	using Sv = ScoreVector<int8_t, SCHAR_MIN>;
	typedef ::DISPATCH_ARCH::SIMD::Vector<int8_t> SeqV;
	constexpr int CHANNELS = ::DISPATCH_ARCH::ScoreTraits<Sv>::CHANNELS;
	assert(subject_count <= CHANNELS);
	
	alignas(CHANNELS) Letter subject_vector[CHANNELS * CHANNELS];
	Sv score, best;
	const Letter* subject_ptr[CHANNELS], * query_end = query + window;
	std::copy(subjects, subjects + subject_count, subject_ptr);

	for (int i = 0; i < window; i += CHANNELS) {
		transpose(subject_ptr, subject_count, subject_vector, SeqV());
		for (size_t j = 0; j < CHANNELS && query < query_end;) {
			SeqV subject_letters(&subject_vector[j * CHANNELS]);
			unsigned query_letter = unsigned(*query);
#ifdef SEQ_MASK
			query_letter &= (unsigned)LETTER_MASK;
#endif
			const Sv match(query_letter, subject_letters);
			score = score + match;
			best = max(best, score);
			++query;
			subject_ptr[j] += CHANNELS;
			++j;
		}
	}

	int8_t best2[CHANNELS];
	best.store(best2);
	const int d = CHANNELS - subject_count;
	for (int i = 0; i < subject_count; ++i)
		out[i] = ScoreTraits<Sv>::int_score(best2[d + i]);
#endif
}

void window_ungapped_best(const Letter* query, const Letter** subjects, int subject_count, int window, int* out) {
#ifdef __SSE4_1__
	if (subject_count < 4) {
#endif
		for (int i = 0; i < subject_count; ++i)
			out[i] = ungapped_window(query, subjects[i], window);
#ifdef __SSE4_1__
	}
#endif
#if ARCH_ID == 2
	else if (subject_count <= 16)
		::DP::ARCH_SSE4_1::window_ungapped(query, subjects, subject_count, window, out);
	else
		::DP::ARCH_AVX2::window_ungapped(query, subjects, subject_count, window, out);
#elif defined(__SSE4_1__)
	window_ungapped(query, subjects, subject_count, window, out);
#endif
}

}}
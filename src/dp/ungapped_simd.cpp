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
#include "../util/simd/transpose16x16.h"

using namespace DISPATCH_ARCH;

namespace DP { namespace DISPATCH_ARCH {

void window_ungapped(const Letter *query, const Letter **subjects, size_t subject_count, int window, int *out) {
#ifdef __SSE4_1__
	typedef score_vector<int8_t> Sv;
	typedef ::DISPATCH_ARCH::SIMD::Vector<int8_t> SeqV;
	assert(window % 16 == 0);

	alignas(16) Letter subject_vector[16 * 16];
	Sv score, best;
	const Letter *subject_ptr[16];
	std::copy(subjects, subjects + subject_count, subject_ptr);

	for (int i = 0; i < window; i += 16) {
		transpose16x16(subject_ptr, subject_count, subject_vector);
		for (int j = 0; j < 16;) {
			SeqV subject_letters(&subject_vector[j * 16]);
			unsigned query_letter = unsigned(*query);
#ifdef SEQ_MASK
			query_letter &= (unsigned)LETTER_MASK;
#endif
			const Sv match(query_letter, subject_letters);
			score = score + match;
			best = max(best, score);
			++query;
			subject_ptr[j] += 16;
			++j;
		}
	}

	int8_t best2[16];
	best.store(best2);
	const size_t d = 16 - subject_count;
	for (size_t i = 0; i < subject_count; ++i)
		out[i] = ScoreTraits<Sv>::int_score(best2[d + i]);
#endif
}

}}
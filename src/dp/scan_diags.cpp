/****
DIAMOND protein aligner
Copyright (C) 2013-2021 Max Planck Society for the Advancement of Science e.V.
						Benjamin Buchfink
						Eberhard Karls Universitaet Tuebingen

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

#include <algorithm>
#include "score_profile.h"
#include "../util/simd.h"
#include "dp.h"
#include "score_vector_int8.h"
#include "../util/simd/dispatch.h"

using namespace DISPATCH_ARCH;

namespace DP { namespace DISPATCH_ARCH {

void scan_diags128(const LongScoreProfile<int8_t>& qp, Sequence s, int d_begin, int j_begin, int j_end, int *out)
{
#ifdef __AVX2__
	using Sv = ScoreVector<int8_t, SCHAR_MIN>;
	const int qlen = (int)qp.length();

	const int j0 = std::max(j_begin, -(d_begin + 128 - 1)),
		i0 = d_begin + j0,
		j1 = std::min(qlen - d_begin, j_end);
	Sv v1, max1, v2, max2, v3, max3, v4, max4;
	for (int i = i0, j = j0; j < j1; ++j, ++i) {
		const int8_t* q = qp.get(s[j], i);
		v1 += Sv(q);
		max1.max(v1);
		q += 32;
		v2 += Sv(q);
		max2.max(v2);
		q += 32;
		v3 += Sv(q);
		max3.max(v3);
		q += 32;
		v4 += Sv(q);
		max4.max(v4);
	}
	int8_t scores[128];
	max1.store(scores);
	max2.store(scores + 32);
	max3.store(scores + 64);
	max4.store(scores + 96);
	for (int i = 0; i < 128; ++i)
		out[i] = ScoreTraits<Sv>::int_score(scores[i]);
#elif defined(__SSE4_1__)
	using Sv = ScoreVector<int8_t, SCHAR_MIN>;
	const int qlen = (int)qp.length();

	const int j0 = std::max(j_begin, -(d_begin + 128 - 1)),
		i0 = d_begin + j0,
		j1 = std::min(qlen - d_begin, j_end);
	Sv v1, max1, v2, max2, v3, max3, v4, max4, v5, max5, v6, max6, v7, max7, v8, max8;
	for (int i = i0, j = j0; j < j1; ++j, ++i) {
		const int8_t* q = qp.get(s[j], i);
		v1 += Sv(q);
		max1.max(v1);
		q += 16;
		v2 += Sv(q);
		max2.max(v2);
		q += 16;
		v3 += Sv(q);
		max3.max(v3);
		q += 16;
		v4 += Sv(q);
		max4.max(v4);
		q += 16;
		v5 += Sv(q);
		max5.max(v5);
		q += 16;
		v6 += Sv(q);
		max6.max(v6);
		q += 16;
		v7 += Sv(q);
		max7.max(v7);
		q += 16;
		v8 += Sv(q);
		max8.max(v8);
	}
	int8_t scores[128];
	max1.store(scores);
	max2.store(scores + 16);
	max3.store(scores + 32);
	max4.store(scores + 48);
	max5.store(scores + 64);
	max6.store(scores + 80);
	max7.store(scores + 96);
	max8.store(scores + 112);
	for (int i = 0; i < 128; ++i)
		out[i] = ScoreTraits<Sv>::int_score(scores[i]);
#else
	const int qlen = (int)qp.length();

	const int j0 = std::max(j_begin, -(d_begin + 128 - 1)),
		i0 = d_begin + j0,
		j1 = std::min(qlen - d_begin, j_end);
	int v[8 * 16], max[8 * 16];
	std::fill(v, v + 128, 0);
	std::fill(max, max + 128, 0);
	for (int i = i0, j = j0; j < j1; ++j, ++i) {
		const int8_t* q = qp.get(s[j], i);
		for (int k = 0; k < 128; ++k) {
			v[k] += q[k];
			v[k] = std::max(v[k], 0);
			max[k] = std::max(max[k], v[k]);
		}
	}
	for (int i = 0; i < 128; ++i)
		out[i] = max[i];
#endif
}

void scan_diags64(const LongScoreProfile<int8_t>& qp, Sequence s, int d_begin, int j_begin, int j_end, int* out)
{
#ifdef __AVX2__
	using Sv = ScoreVector<int8_t, SCHAR_MIN>;
	const int qlen = (int)qp.length();

	const int j0 = std::max(j_begin, -(d_begin + 64 - 1)),
		i0 = d_begin + j0,
		j1 = std::min(qlen - d_begin, j_end);
	Sv v1, max1, v2, max2;
	for (int i = i0, j = j0; j < j1; ++j, ++i) {
		const int8_t* q = qp.get(s[j], i);
		v1 += Sv(q);
		max1.max(v1);
		q += 32;
		v2 += Sv(q);
		max2.max(v2);
	}
	int8_t scores[64];
	max1.store(scores);
	max2.store(scores + 32);
	for (int i = 0; i < 64; ++i)
		out[i] = ScoreTraits<Sv>::int_score(scores[i]);
#elif defined(__SSE4_1__)
	using Sv = ScoreVector<int8_t, SCHAR_MIN>;
	const int qlen = (int)qp.length();

	const int j0 = std::max(j_begin, -(d_begin + 64 - 1)),
		i0 = d_begin + j0,
		j1 = std::min(qlen - d_begin, j_end);
	Sv v1, max1, v2, max2, v3, max3, v4, max4;
	for (int i = i0, j = j0; j < j1; ++j, ++i) {
		const int8_t* q = qp.get(s[j], i);
		v1 += Sv(q);
		max1.max(v1);
		q += 16;
		v2 += Sv(q);
		max2.max(v2);
		q += 16;
		v3 += Sv(q);
		max3.max(v3);
		q += 16;
		v4 += Sv(q);
		max4.max(v4);
	}
	int8_t scores[64];
	max1.store(scores);
	max2.store(scores + 16);
	max3.store(scores + 32);
	max4.store(scores + 48);
	for (int i = 0; i < 64; ++i)
		out[i] = ScoreTraits<Sv>::int_score(scores[i]);
#else
	const int qlen = (int)qp.length();

	const int j0 = std::max(j_begin, -(d_begin + 64 - 1)),
		i0 = d_begin + j0,
		j1 = std::min(qlen - d_begin, j_end);
	int v[64], max[64];
	std::fill(v, v + 64, 0);
	std::fill(max, max + 64, 0);
	for (int i = i0, j = j0; j < j1; ++j, ++i) {
		const int8_t* q = qp.get(s[j], i);
		for (int k = 0; k < 64; ++k) {
			v[k] += q[k];
			v[k] = std::max(v[k], 0);
			max[k] = std::max(max[k], v[k]);
		}
	}
	for (int i = 0; i < 64; ++i)
		out[i] = max[i];
#endif
}

void scan_diags(const LongScoreProfile<int8_t>& qp, Sequence s, int d_begin, int d_end, int j_begin, int j_end, int* out)
{
#ifdef __AVX2__
	using Sv = ScoreVector<int8_t, SCHAR_MIN>;
	const int qlen = (int)qp.length(), band = d_end - d_begin;
	assert(band % 32 == 0);

	const int j0 = std::max(j_begin, -(d_end - 1)),
		i0 = d_begin + j0,
		j1 = std::min(qlen - d_begin, j_end);
	Sv v1, max1, v2, max2;
	for (int i = i0, j = j0; j < j1; ++j, ++i) {
		const int8_t* q = qp.get(s[j], i);
		v1 += Sv(q);
		max1.max(v1);
		q += 32;
		v2 += Sv(q);
		max2.max(v2);
	}
	int8_t scores[64];
	max1.store(scores);
	max2.store(scores + 32);
	for (int i = 0; i < 64; ++i)
		out[i] = ScoreTraits<Sv>::int_score(scores[i]);
#elif defined(__SSE4_1__)
	using Sv = ScoreVector<int8_t, SCHAR_MIN>;
	const int qlen = (int)qp.length();

	const int j0 = std::max(j_begin, -(d_begin + 64 - 1)),
		i0 = d_begin + j0,
		j1 = std::min(qlen - d_begin, j_end);
	Sv v1, max1, v2, max2, v3, max3, v4, max4;
	for (int i = i0, j = j0; j < j1; ++j, ++i) {
		const int8_t* q = qp.get(s[j], i);
		v1 += Sv(q);
		max1.max(v1);
		q += 16;
		v2 += Sv(q);
		max2.max(v2);
		q += 16;
		v3 += Sv(q);
		max3.max(v3);
		q += 16;
		v4 += Sv(q);
		max4.max(v4);
	}
	int8_t scores[64];
	max1.store(scores);
	max2.store(scores + 16);
	max3.store(scores + 32);
	max4.store(scores + 48);
	for (int i = 0; i < 64; ++i)
		out[i] = ScoreTraits<Sv>::int_score(scores[i]);
#else
	const int qlen = (int)qp.length();

	const int j0 = std::max(j_begin, -(d_begin + 64 - 1)),
		i0 = d_begin + j0,
		j1 = std::min(qlen - d_begin, j_end);
	int v[64], max[64];
	std::fill(v, v + 64, 0);
	std::fill(max, max + 64, 0);
	for (int i = i0, j = j0; j < j1; ++j, ++i) {
		const int8_t* q = qp.get(s[j], i);
		for (int k = 0; k < 64; ++k) {
			v[k] += q[k];
			v[k] = std::max(v[k], 0);
			max[k] = std::max(max[k], v[k]);
		}
	}
	for (int i = 0; i < 64; ++i)
		out[i] = max[i];
#endif
}

int diag_alignment(const int* s, int count) {
	int best = 0, best_gap = -score_matrix.gap_open(), d = -1;
	for (int i = 0; i < count; ++i) {
		if (s[i] < config.gapped_filter_diag_score)
			continue;
		const int gap_score = -score_matrix.gap_extend() * (i - d) + best_gap;
		int n = s[i];
		if (gap_score + s[i] > best) {
			best = n = gap_score + s[i];
		}
		if (s[i] > best) {
			best = n = s[i];
		}
		const int open_score = -score_matrix.gap_open() + n;
		if (open_score > gap_score) {
			best_gap = open_score;
			d = i;
		}
	}
	return best;
}

}

DISPATCH_6V(scan_diags128, const LongScoreProfile<int8_t>&, qp, Sequence, s, int, d_begin, int, j_begin, int, j_end, int*, out)
DISPATCH_6V(scan_diags64, const LongScoreProfile<int8_t>&, qp, Sequence, s, int, d_begin, int, j_begin, int, j_end, int*, out)
DISPATCH_7V(scan_diags, const LongScoreProfile<int8_t>&, qp, Sequence, s, int, d_begin, int, d_end, int, j_begin, int, j_end, int*, out)
DISPATCH_2(int, diag_alignment, const int*, s, int, count)

}
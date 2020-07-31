/****
DIAMOND protein aligner
Copyright (C) 2016-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
						
Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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

#pragma once
#include <assert.h>
#include "../score_vector.h"
#include "../score_vector_int8.h"
#include "../score_vector_int16.h"
#include "../../basic/value.h"
#include "../util/simd/vector.h"

static inline uint8_t cmp_mask(int x, int y) {
	return x == y;
}

static inline int blend(int v, int w, int mask) {
	return mask ? w : v;
}

template<typename _sv>
struct TraceStat {
	_sv length;
	_sv gapopen;
	_sv qstart;
	_sv sstart;
	_sv ident;
	_sv mismatch;
};

template<typename _sv>
struct RowCounter {
	RowCounter(int offset):
		i(::DISPATCH_ARCH::ScoreTraits<_sv>::zero_score() + typename ::DISPATCH_ARCH::ScoreTraits<_sv>::Score(offset)),
		i_max()
	{
	}
	_sv i;
	_sv i_max;
};

template<typename _sv>
static inline _sv swipe_cell_update(const _sv &diagonal_cell,
	const _sv &scores,
	void*,
	const _sv &gap_extension,
	const _sv &gap_open,
	_sv &horizontal_gap,
	_sv &vertical_gap,
	_sv &best,
	void*,
	void*,
	void*,
	void*,
	const RowCounter<_sv>&)
{
	using std::max;
	_sv current_cell = max(diagonal_cell + scores, vertical_gap);
	current_cell = max(current_cell, horizontal_gap);
	::DISPATCH_ARCH::ScoreTraits<_sv>::saturate(current_cell);
	best = max(best, current_cell);
	vertical_gap -= gap_extension;
	horizontal_gap -= gap_extension;
	const _sv open = current_cell - gap_open;
	vertical_gap = max(vertical_gap, open);
	horizontal_gap = max(horizontal_gap, open);
	return current_cell;
}


template<typename _sv>
static inline _sv swipe_cell_update(const _sv &diagonal_cell,
	const _sv &scores,
	void*,
	const _sv &gap_extension,
	const _sv &gap_open,
	_sv &horizontal_gap,
	_sv &vertical_gap,
	_sv &best,
	void*,
	void*,
	void*,
	typename ::DISPATCH_ARCH::ScoreTraits<_sv>::TraceMask *trace_mask,
	RowCounter<_sv>& row_counter)
{
	typedef typename ::DISPATCH_ARCH::ScoreTraits<_sv>::TraceMask TraceMask;
	using std::max;
	_sv current_cell = max(diagonal_cell + scores, vertical_gap);
	current_cell = max(current_cell, horizontal_gap);
	::DISPATCH_ARCH::ScoreTraits<_sv>::saturate(current_cell);

	trace_mask->gap = TraceMask::make(cmp_mask(current_cell, vertical_gap), cmp_mask(current_cell, horizontal_gap));

	best = max(best, current_cell);
	row_counter.i_max = blend(row_counter.i_max, row_counter.i, best == current_cell);
	row_counter.i += _sv(typename ::DISPATCH_ARCH::ScoreTraits<_sv>::Score(1));

	vertical_gap -= gap_extension;
	horizontal_gap -= gap_extension;
	const _sv open = current_cell - gap_open;
	vertical_gap = max(vertical_gap, open);
	horizontal_gap = max(horizontal_gap, open);

	trace_mask->open = TraceMask::make(cmp_mask(vertical_gap, open), cmp_mask(horizontal_gap, open));

	return current_cell;
}

template<typename _sv>
static inline _sv swipe_cell_update(const _sv &diagonal_cell,
	const _sv &scores,
	void*,
	const _sv &gap_extension,
	const _sv &gap_open,
	_sv &horizontal_gap,
	_sv &vertical_gap,
	_sv &best,
	TraceStat<_sv> &trace_stat_diag,
	TraceStat<_sv> &trace_stat_vertical,
	TraceStat<_sv> &trace_stat_horizontal,
	void*,
	const RowCounter<_sv>&)
{
	using std::max;
	_sv current_cell = max(diagonal_cell + scores, vertical_gap);
	current_cell = max(current_cell, horizontal_gap);
	::DISPATCH_ARCH::ScoreTraits<_sv>::saturate(current_cell);
	best = max(best, current_cell);
	vertical_gap -= gap_extension;
	horizontal_gap -= gap_extension;
	const _sv open = current_cell - gap_open;
	vertical_gap = max(vertical_gap, open);
	horizontal_gap = max(horizontal_gap, open);
	return current_cell;
}

template<typename _sv>
static inline _sv swipe_cell_update(const _sv& diagonal_cell,
	const _sv& scores,
	const _sv& query_bias,
	const _sv& gap_extension,
	const _sv& gap_open,
	_sv& horizontal_gap,
	_sv& vertical_gap,
	_sv& best,
	void*,
	void*,
	void*,
	void*,
	const RowCounter<_sv>&)
{
	using std::max;
	_sv current_cell = diagonal_cell + (scores + query_bias);
	current_cell = max(max(current_cell, vertical_gap), horizontal_gap);
	::DISPATCH_ARCH::ScoreTraits<_sv>::saturate(current_cell);
	best = max(best, current_cell);
	vertical_gap -= gap_extension;
	horizontal_gap -= gap_extension;
	const _sv open = current_cell - gap_open;
	vertical_gap = max(vertical_gap, open);
	horizontal_gap = max(horizontal_gap, open);
	return current_cell;
}

template<typename _sv>
static inline _sv swipe_cell_update(const _sv& diagonal_cell,
	const _sv& scores,
	const _sv& query_bias,
	const _sv& gap_extension,
	const _sv& gap_open,
	_sv& horizontal_gap,
	_sv& vertical_gap,
	_sv& best,
	TraceStat<_sv> &trace_stat_diag,
	TraceStat<_sv> &trace_stat_vertical,
	TraceStat<_sv> &trace_stat_horizontal,
	void*,
	const RowCounter<_sv>&)
{
	typedef typename ::DISPATCH_ARCH::ScoreTraits<_sv>::Score Score;
	using std::max;
	_sv current_cell = diagonal_cell + (scores + query_bias);
	current_cell = max(max(current_cell, vertical_gap), horizontal_gap);
	::DISPATCH_ARCH::ScoreTraits<_sv>::saturate(current_cell);

	const _sv one = _sv(Score(1)), zero = _sv(), zero2 = _sv(Score(0));
	const _sv vgap_mask = current_cell == vertical_gap, hgap_mask = current_cell == horizontal_gap, zero_mask = current_cell == zero;

	best = max(best, current_cell);
	vertical_gap -= gap_extension;
	horizontal_gap -= gap_extension;
	const _sv open = current_cell - gap_open;
	vertical_gap = max(vertical_gap, open);
	horizontal_gap = max(horizontal_gap, open);

	/*trace_stat_vertical.length += one;
	trace_stat_horizontal.length += one;
	trace_stat_diag.length += one;
	trace_stat_diag.length = blend(trace_stat_diag.length, trace_stat_vertical.length, vgap_mask);
	trace_stat_diag.length = blend(trace_stat_diag.length, trace_stat_horizontal.length, hgap_mask);
	trace_stat_diag.length = blend(trace_stat_diag.length, zero2, zero_mask);*/
	
	return current_cell;
}

template<typename _sv>
static inline _sv swipe_cell_update(const _sv& diagonal_cell,
	const _sv& scores,
	const _sv& query_bias,
	const _sv& gap_extension,
	const _sv& gap_open,
	_sv& horizontal_gap,
	_sv& vertical_gap,
	_sv& best,
	void*,
	void*,
	void*,
	typename ::DISPATCH_ARCH::ScoreTraits<_sv>::TraceMask *trace_mask,
	RowCounter<_sv> &row_counter)
{
	typedef typename ::DISPATCH_ARCH::ScoreTraits<_sv>::TraceMask TraceMask;
	using std::max;
	_sv current_cell = diagonal_cell + (scores + query_bias);
	current_cell = max(max(current_cell, vertical_gap), horizontal_gap);
	::DISPATCH_ARCH::ScoreTraits<_sv>::saturate(current_cell);

	trace_mask->gap = TraceMask::make(cmp_mask(current_cell, vertical_gap), cmp_mask(current_cell, horizontal_gap));

	best = max(best, current_cell);
	row_counter.i_max = blend(row_counter.i_max, row_counter.i, best == current_cell);
	row_counter.i += _sv(typename ::DISPATCH_ARCH::ScoreTraits<_sv>::Score(1));

	vertical_gap -= gap_extension;
	horizontal_gap -= gap_extension;
	const _sv open = current_cell - gap_open;
	vertical_gap = max(vertical_gap, open);
	horizontal_gap = max(horizontal_gap, open);

	trace_mask->open = TraceMask::make(cmp_mask(vertical_gap, open), cmp_mask(horizontal_gap, open));

	return current_cell;
}

template<typename _sv>
static inline _sv cell_update(const _sv &diagonal_cell,
	const _sv &scores,
	const _sv &gap_extension,
	const _sv &gap_open,
	_sv &horizontal_gap,
	_sv &vertical_gap,
	_sv &best,
	const _sv &vbias)
{
	_sv current_cell = diagonal_cell + scores;
	current_cell -= vbias;
	current_cell.max(vertical_gap).max(horizontal_gap);
	best.max(current_cell);
	vertical_gap -= gap_extension;
	horizontal_gap -= gap_extension;
	const _sv open = current_cell - gap_open;
	vertical_gap.max(open);
	horizontal_gap.max(open);
	return current_cell;
}

template<typename _score>
static inline DISPATCH_ARCH::score_vector<_score> cell_update(const DISPATCH_ARCH::score_vector<_score> &diagonal_cell,
	const DISPATCH_ARCH::score_vector<_score> &scores,
	const DISPATCH_ARCH::score_vector<_score> &gap_extension,
	const DISPATCH_ARCH::score_vector<_score> &gap_open,
	DISPATCH_ARCH::score_vector<_score> &horizontal_gap,
	DISPATCH_ARCH::score_vector<_score> &vertical_gap,
	DISPATCH_ARCH::score_vector<_score> &best)
{
	using namespace ::DISPATCH_ARCH;
	score_vector<_score> current_cell = diagonal_cell + scores;
	current_cell.max(vertical_gap).max(horizontal_gap);
	best.max(current_cell);
	vertical_gap -= gap_extension;
	horizontal_gap -= gap_extension;
	const score_vector<_score> open = current_cell - gap_open;
	vertical_gap.max(open);
	horizontal_gap.max(open);
	return current_cell;
}

template<typename _sv>
static inline _sv cell_update(const _sv &diagonal_cell,
	const _sv &shift_cell0,
	const _sv &shift_cell1,
	const _sv &scores,
	const _sv &gap_extension,
	const _sv &gap_open,
	const _sv &frame_shift,
	_sv &horizontal_gap,
	_sv &vertical_gap,
	_sv &best)
{
	using std::max;
	_sv current_cell = diagonal_cell + scores;
	const _sv f = scores - frame_shift;
	current_cell = max(current_cell, shift_cell0 + f);
	current_cell = max(current_cell, shift_cell1 + f);
	current_cell = max(max(current_cell, vertical_gap), horizontal_gap);
	::DISPATCH_ARCH::ScoreTraits<_sv>::saturate(current_cell);
	best = max(best, current_cell);
	vertical_gap -= gap_extension;
	horizontal_gap -= gap_extension;
	const _sv open = current_cell - gap_open;
	vertical_gap = max(vertical_gap, open);
	horizontal_gap = max(horizontal_gap, open);
	return current_cell;
}

namespace DISPATCH_ARCH {

template<typename _sv>
struct SwipeProfile
{

#ifdef __SSSE3__
	inline void set(typename ScoreTraits<_sv>::Vector seq)
	{
		assert(sizeof(data_) / sizeof(_sv) >= value_traits.alphabet_size);
		for (unsigned j = 0; j < AMINO_ACID_COUNT; ++j)
			data_[j] = _sv(j, seq);
	}
#else
	inline void set(uint64_t seq)
	{
		assert(sizeof(data_) / sizeof(_sv) >= value_traits.alphabet_size);
		for (unsigned j = 0; j < AMINO_ACID_COUNT; ++j)
			data_[j] = _sv(j, seq);
	}
#endif
	inline const _sv& get(Letter i) const
	{
		return data_[(int)i];
	}
	_sv data_[AMINO_ACID_COUNT];
};

template<>
struct SwipeProfile<int32_t>
{
#ifdef __AVX2__
	void set(const __m256i& seq)
	{
		int16_t s[32];
		_mm256_storeu_si256((__m256i*)s, seq);
		row = score_matrix.row((char)s[0]);
	}
#endif
#ifdef __SSE2__
	void set(const __m128i& seq)
	{
		int16_t s[8];
		_mm_storeu_si128((__m128i*)s, seq);
		row = score_matrix.row((char)s[0]);
	}
#endif
	void set(uint64_t seq)
	{
		row = score_matrix.row((char)seq);
	}
	int32_t get(char i) const
	{
		return row[(int)i];
	}
	const int32_t *row;
};

}

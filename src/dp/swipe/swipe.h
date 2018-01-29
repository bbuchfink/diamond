/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#ifndef SWIPE_H_
#define SWIPE_H_

#include "../score_vector.h"
#include "../../basic/value.h"

template<typename _sv>
inline _sv cell_update(const _sv &diagonal_cell,
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
inline score_vector<_score> cell_update(const score_vector<_score> &diagonal_cell,
	const score_vector<_score> &scores,
	const score_vector<_score> &gap_extension,
	const score_vector<_score> &gap_open,
	score_vector<_score> &horizontal_gap,
	score_vector<_score> &vertical_gap,
	score_vector<_score> &best)
{
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
inline _sv cell_update(const _sv &diagonal_cell,
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
	ScoreTraits<_sv>::saturate(current_cell);
	best = max(best, current_cell);
	vertical_gap -= gap_extension;
	horizontal_gap -= gap_extension;
	const _sv open = current_cell - gap_open;
	vertical_gap = max(vertical_gap, open);
	horizontal_gap = max(horizontal_gap, open);
	return current_cell;
}

template<typename _sv>
struct SwipeProfile
{
#ifdef __SSSE3__
	inline void set(const __m128i &seq)
	{
		assert(sizeof(data_) / sizeof(_sv) >= value_traits.alphabet_size);
		_sv bias(score_matrix.bias());
		for (unsigned j = 0; j < AMINO_ACID_COUNT; ++j)
			data_[j] = _sv(j, seq, bias);
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

#endif
/****
DIAMOND protein aligner
Copyright (C) 2016-2021 Max Planck Society for the Advancement of Science e.V.
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
#include <limits>
#include <vector>
#include "../score_vector.h"
#include "../score_vector_int8.h"
#include "../score_vector_int16.h"
#include "../../basic/value.h"
#include "../../util/simd/vector.h"
#include "../../util/system.h"
#include "../../util/memory/alignment.h"
#include "../../util/simd/transpose.h"

namespace DP {
	struct NoCBS;
}

template<typename _sv, typename _cbs>
struct CBSBuffer {
	CBSBuffer(const DP::NoCBS&, int, uint32_t) {}
	void* operator()(int i) const {
		return nullptr;
	}
};

template<typename _sv>
struct CBSBuffer<_sv, const int8_t*> {
	CBSBuffer(const int8_t* v, int l, uint32_t channel_mask) {
		typedef typename ::DISPATCH_ARCH::ScoreTraits<_sv>::Score Score;
		data.reserve(l);
		for (int i = 0; i < l; ++i)
			data.push_back(blend_sv<_sv>(Score(v[i]), (Score)0, channel_mask));
	}
	_sv operator()(int i) const {
		return data[i];
	}
	std::vector<_sv, Util::Memory::AlignmentAllocator<_sv, 32>> data;
};

template<typename _sv>
FORCE_INLINE _sv cell_update(const _sv &diagonal_cell,
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
	saturate(current_cell);
	best = max(best, current_cell);
	vertical_gap -= gap_extension;
	horizontal_gap -= gap_extension;
	const _sv open = current_cell - gap_open;
	vertical_gap = max(vertical_gap, open);
	horizontal_gap = max(horizontal_gap, open);
	return current_cell;
}

namespace DISPATCH_ARCH {

template<typename Sv>
struct SwipeProfile
{

#ifdef __SSSE3__
	inline void set(typename ScoreTraits<Sv>::Vector seq)
	{
		assert(sizeof(data_) / sizeof(Sv) >= value_traits.alphabet_size);
		for (unsigned j = 0; j < AMINO_ACID_COUNT; ++j)
			data_[j] = Sv(j, seq);
	}
#endif

	inline const Sv& get(Letter i) const
	{
		return data_[(int)i];
	}

	void set(const int8_t** target_scores) {
#if ARCH_ID == 2
		transpose(target_scores, 32, (int8_t*)data_, __m256i());
		for (size_t i = 0; i < AMINO_ACID_COUNT; ++i)
			data_[i].expand_from_8bit();
#elif defined(__SSE2__)
		transpose(target_scores, 16, (int8_t*)data_, __m128i());
		for (int i = 0; i < 16; ++i)
			target_scores[i] += 16;
		transpose(target_scores, 16, (int8_t*)(&data_[16]), __m128i());
		for (size_t i = 0; i < AMINO_ACID_COUNT; ++i)
			data_[i].expand_from_8bit();
#else
		for (int i = 0; i < AMINO_ACID_COUNT; ++i)
			data_[i] = target_scores[0][i];
#endif
	}

	void set(const int32_t** target_scores) {
		typename ScoreTraits<Sv>::Score s[ScoreTraits<Sv>::CHANNELS];
		for (size_t i = 0; i < AMINO_ACID_COUNT; ++i) {
			for (size_t j = 0; j < ScoreTraits<Sv>::CHANNELS; ++j)
				s[j] = target_scores[j][i];
			data_[i] = load_sv<Sv>(s);
		}
	}

	//_sv data_[AMINO_ACID_COUNT];
	Sv data_[32];

};

template<>
struct SwipeProfile<int32_t>
{
#ifdef __AVX2__
	void set(const __m256i& seq)
	{
		int16_t s[32];
		_mm256_storeu_si256((__m256i*)s, seq);
		const int* row = score_matrix.row((char)s[0]);
		std::copy(row, row + 32, this->row);
	}
#endif
#ifdef __SSE2__
	void set(const __m128i& seq)
	{
		int16_t s[8];
		_mm_storeu_si128((__m128i*)s, seq);
		const int* row = score_matrix.row((char)s[0]);
		std::copy(row, row + 32, this->row);
	}
#endif
	void set(uint64_t seq)
	{
		const int* row = score_matrix.row((char)seq);
		std::copy(row, row + 32, this->row);
	}
	void set(const int8_t** target_scores) {
		for (int i = 0; i < 32; ++i)
			row[i] = target_scores[0][i];
	}
	void set(const int32_t** target_scores) {
		for (int i = 0; i < 32; ++i)
			row[i] = target_scores[0][i];
	}
	int32_t get(char i) const
	{
		return row[(int)i];
	}
	int32_t row[32];
};

}

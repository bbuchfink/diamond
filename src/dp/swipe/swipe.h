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

#pragma once
#include <assert.h>
#include <vector>
#include "../score_vector.h"
#include "basic/value.h"
#include "util/system.h"
#include "util/memory/alignment.h"
#include "util/simd/transpose.h"
#include "stats/score_matrix.h"

namespace DP {
	struct NoCBS;
}

template<typename Sv, typename Cbs>
struct CBSBuffer {
	CBSBuffer(const DP::NoCBS&, int, uint32_t) {}
	void* operator()(int i) const {
		return nullptr;
	}
};

template<typename Sv>
struct CBSBuffer<Sv, const int8_t*> {
	CBSBuffer(const int8_t* v, int l, uint32_t channel_mask) {
		typedef typename ::DISPATCH_ARCH::ScoreTraits<Sv>::Score Score;
		data.reserve(l);
		for (int i = 0; i < l; ++i)
			data.push_back(blend_sv<Sv>(Score(v[i]), (Score)0, channel_mask));
	}
	Sv operator()(int i) const {
		return data[i];
	}
	std::vector<Sv, Util::Memory::AlignmentAllocator<Sv, 32>> data;
};

template<typename Sv>
FORCE_INLINE Sv cell_update(const Sv &diagonal_cell,
	const Sv &shift_cell0,
	const Sv &shift_cell1,
	const Sv &scores,
	const Sv &gap_extension,
	const Sv &gap_open,
	const Sv &frame_shift,
	Sv &horizontal_gap,
	Sv &vertical_gap,
	Sv &best)
{
	using std::max;
	Sv current_cell = diagonal_cell + scores;
	const Sv f = scores - frame_shift;
	current_cell = max(current_cell, shift_cell0 + f);
	current_cell = max(current_cell, shift_cell1 + f);
	current_cell = max(max(current_cell, vertical_gap), horizontal_gap);
	saturate(current_cell);
	best = max(best, current_cell);
	vertical_gap -= gap_extension;
	horizontal_gap -= gap_extension;
	const Sv open = current_cell - gap_open;
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
#elif defined(__ARM_NEON)
		transpose(target_scores, 16, (int8_t*)data_, int8x16_t());
		for (int i = 0; i < 16; ++i)
			target_scores[i] += 16;
		transpose(target_scores, 16, (int8_t*)(&data_[16]), int8x16_t());
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
#ifdef __ARM_NEON
	void set(const int16x8_t& seq)
	{
		int16_t s[8];
		vst1q_s16(s, seq);
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
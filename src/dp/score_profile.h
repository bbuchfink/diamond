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
#include <vector>
#include "../basic/sequence.h"
#include "../basic/value.h"
#include "../stats/hauser_correction.h"

template<typename Score = int8_t>
struct LongScoreProfile
{
	enum { DEFAULT_PADDING = 128 };
	LongScoreProfile(int64_t padding = DEFAULT_PADDING):
		padding(std::max(padding, (int64_t)DEFAULT_PADDING))
	{}
	size_t length() const
	{
		return data[0].size() - 2 * padding;
	}
	const Score* get(Letter l, int i) const
	{
		return &data[(int)l][i + padding];
	}
	std::vector<const Score*> pointers(int offset) const {
		std::vector<const Score*> v;
		v.reserve(AMINO_ACID_COUNT);
		for (size_t i = 0; i < AMINO_ACID_COUNT; ++i)
			v.push_back(get(Letter(i), offset));
		return v;
	}
	LongScoreProfile reverse() const {
		LongScoreProfile r(*this);
		for (size_t i = 0; i < AMINO_ACID_COUNT; ++i)
			std::reverse(r.data[i].begin(), r.data[i].end());
		return r;
	}
	std::vector<Score> data[AMINO_ACID_COUNT];
	int64_t padding;
};

namespace DP {

LongScoreProfile<int8_t> make_profile8(Sequence seq, const int8_t* cbs, int64_t padding);
LongScoreProfile<int16_t> make_profile16(Sequence seq, const int8_t* cbs, int64_t padding);

}
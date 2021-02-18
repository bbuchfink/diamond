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
#include "score_vector.h"
#include "../basic/value.h"
#include "../stats/hauser_correction.h"

struct LongScoreProfile
{
	LongScoreProfile()
	{}
	LongScoreProfile(Sequence seq, const Bias_correction &cbs)
	{
		for (unsigned l = 0; l < AMINO_ACID_COUNT; ++l) {
			const int8_t* scores = &score_matrix.matrix8()[l << 5];
			data[l].reserve(seq.length() + 2 * padding);
			data[l].insert(data[l].end(), padding, 0);
			for (unsigned i = 0; i < seq.length(); ++i)
				data[l].push_back(scores[(int)seq[i]] + cbs.int8[i]);
			data[l].insert(data[l].end(), padding, 0);
		}
	}
	LongScoreProfile(Sequence seq)
	{
		set(seq);
	}
	void set(Sequence seq) {
		for (unsigned l = 0; l < AMINO_ACID_COUNT; ++l) {
			const int8_t* scores = &score_matrix.matrix8()[l << 5];
			data[l].clear();
			data[l].reserve(seq.length() + 2 * padding);
			data[l].insert(data[l].end(), padding, 0);
			for (unsigned i = 0; i < seq.length(); ++i)
				data[l].push_back(scores[(int)seq[i]]);
			data[l].insert(data[l].end(), padding, 0);
		}
	}
	size_t length() const
	{
		return data[0].size() - 2 * padding;
	}
	const int8_t* get(Letter l, int i) const
	{
		return &data[(int)l][i + padding];
	}
	std::vector<int8_t> data[AMINO_ACID_COUNT];
	enum { padding = 128 };
};

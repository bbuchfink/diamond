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
#include <vector>
#include <stdint.h>
#include "../basic/sequence.h"
#include "../util/geo/diagonal_segment.h"
#include "../basic/match.h"

struct No_score_correction
{
	void operator()(int &score, int i, int query_anchor, int mult) const
	{}
};

struct Bias_correction : public std::vector<float>
{
	Bias_correction() {}
	Bias_correction(const Sequence &seq);
	void operator()(float &score, int i, int query_anchor, int mult) const
	{
		score += (*this)[query_anchor + i * mult];
	}
	int operator()(const Hsp &hsp) const;
	int operator()(const DiagonalSegment &d) const;
	static std::vector<int8_t> reverse(const int8_t* p, const size_t len);
	std::vector<int8_t> int8;
};

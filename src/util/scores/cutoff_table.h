/****
DIAMOND protein aligner
Copyright (C) 2020 Max Planck Society for the Advancement of Science e.V.

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
#include "../../stats/score_matrix.h"
#include "../intrin.h"

namespace Util { namespace Scores {

struct CutoffTable {

	CutoffTable() {}
	
	CutoffTable(double evalue) {
		for (int b = 1; b <= MAX_BITS; ++b) {
			data_[b] = score_matrix.rawscore(score_matrix.bitscore_norm(evalue, 1 << (b - 1)));
		}
	}

	int operator()(int query_len) const {
		const int b = 32 - clz((uint32_t)query_len);
		return data_[b];
	}

private:

	enum { MAX_BITS = 31 };

	int data_[MAX_BITS + 1];

};

struct CutoffTable2D {

	CutoffTable2D() {}

	CutoffTable2D(double evalue) {
		for (int b1 = 1; b1 <= MAX_BITS; ++b1)
			for (int b2 = 1; b2 <= MAX_BITS; ++b2) {
				data_[b1][b2] = calc_min_score(1 << (b1 - 1), 1 << (b2 - 1), evalue);
			}
	}

	int operator()(int query_len, int target_len) const {
		const int b1 = 32 - clz((uint32_t)query_len), b2 = 32 - clz((uint32_t)target_len);
		return data_[b1][b2];
	}

private:

	int calc_min_score(unsigned qlen, unsigned slen, double evalue) {
		for (int i = 10; i < 1000; ++i)
			if (score_matrix.evalue_norm(i, qlen, slen) <= evalue)
				return i;
		return 1000;
	}

	enum { MAX_BITS = 31 };

	int data_[MAX_BITS + 1][MAX_BITS + 1];

};

}}
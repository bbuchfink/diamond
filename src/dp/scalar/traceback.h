/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

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
#include "../basic/match.h"
#include "../stats/score_matrix.h"

template<typename _matrix>
bool have_vgap(const _matrix& dp,
	int i,
	int j,
	int gap_open,
	int gap_extend,
	int& l)
{
	int score = dp(i, j);
	l = 1;
	--i;
	while (i > 0) {
		if (score == dp(i, j) - gap_open - l * gap_extend)
			return true;
		--i;
		++l;
	}
	return false;
}

template<typename _matrix>
bool have_hgap(const _matrix& dp,
	int i,
	int j,
	int gap_open,
	int gap_extend,
	int& l)
{
	int score = dp(i, j);
	l = 1;
	--j;
	while (j > 0) {
		if (score == dp(i, j) - gap_open - l * gap_extend)
			return true;
		--j;
		++l;
	}
	return false;
}

template<typename _matrix>
int have_diag(const _matrix& dp,
	int i,
	int j,
	const Sequence& query,
	const Sequence& subject,
	bool log)
{
	int l = 0;
	while (i > 0 && j > 0) {
		const int match_score = score_matrix(query[i - 1], subject[j - 1]);

		if (dp(i, j) == match_score + dp(i - 1, j - 1)) {
			/*if (log)
				printf("i=%i j=%i score=%i subject=%c query=%c\n", i, j, dp(i, j), value_traits.alphabet[(int)subject[j - 1]], value_traits.alphabet[(int)query[i - 1]]);*/
			++l;
			--i;
			--j;
		}
		else
			break;
	}
	return l;
}

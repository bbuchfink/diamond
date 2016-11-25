/****
Copyright (c) 2016, Benjamin Buchfink
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/

#ifndef TRACEBACK_H_
#define TRACEBACK_H_

#include "../basic/match.h"

template<typename _matrix>
bool have_vgap(const _matrix &dp,
	int i,
	int j,
	int gap_open,
	int gap_extend,
	int &l)
{
	int score = dp(i, j);
	l = 1;
	--i;
	while (i > 0) {
		if (score == dp(i, j) - gap_open - l*gap_extend)
			return true;
		--i;
		++l;
	}
	return false;
}

template<typename _matrix>
bool have_hgap(const _matrix &dp,
	int i,
	int j,
	int gap_open,
	int gap_extend,
	int &l)
{
	int score = dp(i, j);
	l = 1;
	--j;
	while (j > 0) {
		if (score == dp(i, j) - gap_open - l*gap_extend)
			return true;
		--j;
		++l;
	}
	return false;
}

template<typename _matrix>
int have_diag(const _matrix &dp,
	int i,
	int j,
	const sequence &query,
	const sequence &subject,
	bool log)
{
	int l = 0;
	while (i > 0 && j > 0) {
		const int match_score = score_matrix(query[i - 1], subject[j - 1]);

		if (dp(i, j) == match_score + dp(i - 1, j - 1)) {
			if (log)
				printf("i=%i j=%i score=%i subject=%c query=%c\n", i, j, dp(i, j), value_traits.alphabet[(int)subject[j - 1]], value_traits.alphabet[(int)query[i - 1]]);
			++l;
			--i;
			--j;
		}
		else
			break;
	}
	return l;
}

#endif
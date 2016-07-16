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
	--j;
	while (dp.in_band(i, j)) {
		if (score == dp(i, j) - gap_open - (l - 1)*gap_extend)
			return true;
		--j;
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
	--i;
	while (dp.in_band(i, j)) {
		if (score == dp(i, j) - gap_open - (l - 1)*gap_extend)
			return true;
		--i;
		++l;
	}
	return false;
}

template<typename _dir, typename _matrix>
Hsp_data traceback(const Letter *query,
	const Letter *subject,
	const _matrix &dp,
	int gap_open,
	int gap_extend,
	int subject_pos,
	int query_pos,
	int score,
	vector<char> &transcript_buf)
{
	if (i == -1)
		return local_match(0);

	Hsp_data l;
	l.score_ = score;
	Edit_transcript transcript(transcript_buf);

	int gap_len, i = subject_pos, j = query_pos;

	while (dp(i,j) != 0) {
		const Letter lq = get_dir(query, j, _dir()), ls = mask_critical(get_dir(subject, i, _dir()));
		const int match_score = score_matrix(lq, ls);
		//printf("i=%i j=%i score=%i subject=%c query=%c\n",i,j,dp(i, j),Value_traits<_val>::ALPHABET[ls],Value_traits<_val>::ALPHABET[lq]);

		if (dp(i, j) == match_score + dp(i - 1, j - 1)) {
			if (lq == ls)
				++l.identities_;
			else
				++l.mismatches_;
			--i;
			--j;
			++l.len_;
			transcript_buf.push_back(op_match);
		}
		else if (have_hgap(dp, i, j, gap_open, gap_extend, gap_len)) {
			++l.gap_openings_;
			l.len_ += gap_len;
			i -= gap_len;
			transcript_buf.insert(transcript_buf.end(), gap_len, op_deletion);
		}
		else if (have_vgap(dp, i, j, gap_open, gap_extend, gap_len)) {
			++l.gap_openings_;
			l.len_ += gap_len;
			j -= gap_len;
			transcript_buf.insert(transcript_buf.end(), gap_len, op_insertion);
		} else
			throw std::runtime_error("Traceback error.");
	}

	l.transcript_right_ = transcript.set_end(transcript_buf);
	return l;

}

#endif
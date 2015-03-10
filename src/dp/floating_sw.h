/****
Copyright (c) 2014, University of Tuebingen
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****
Author: Benjamin Buchfink
****/

#ifndef FLOATING_SW_H_
#define FLOATING_SW_H_

#include "../basic/match.h"
#include "scalar_dp_matrix.h"
#include "../util/direction.h"
#include "scalar_traceback.h"

template<typename _val, typename _dir, typename _score, typename _traceback>
local_match<_val> floating_sw_dir(const _val *query, const _val *subject, int band, _score xdrop, _score gap_open, _score gap_extend)
{
	using std::max;

	_score max_score = 0, column_max = 0;
	int j = 0, i_max = -1, j_best = -1, i_best = -1;
	Scalar_dp_matrix<_score,_traceback> mtx (band);
	const _val *x = query, *y = subject;

	while(*y != String_set<_val>::PADDING_CHAR && max_score - column_max < xdrop) {
		typename Scalar_dp_matrix<_score,_traceback>::Column_iterator it = mtx.column(j, i_max);
		if(get_dir(x, it.row(), _dir()) == String_set<_val>::PADDING_CHAR)
			break;
		_score vgap = Scalar_dp_matrix<_score,_traceback>::NEG_MIN;
		if(get_dir(x, i_max+1, _dir()) == String_set<_val>::PADDING_CHAR) {
			column_max = std::numeric_limits<_score>::min();
		} else {
			++i_max;
			column_max += score_matrix::get().letter_score(mask_critical(*y), get_dir(x, i_max, _dir()));
		}

		for(; it.valid() && get_dir(x, it.row(), _dir()) != String_set<_val>::PADDING_CHAR; ++it) {
			const _score match_score = score_matrix::get().letter_score(mask_critical(*y), get_dir(x, it.row(), _dir()));
			const _score s = max(max(it.diag() + match_score, vgap), it.hgap_in());
			if(s > column_max) {
				column_max = s;
				i_max = it.row();
			}
			const _score open = s - gap_open;
			vgap = max(vgap - gap_extend, open);
			it.hgap_out() = max(it.hgap_in() - gap_extend, open);
			it.score() = s;
		}

		if(column_max > max_score) {
			max_score = column_max;
			j_best = j;
			i_best = i_max;
		}
		y = inc_dir(y, _dir());
		++j;
	}

	return traceback<_val,_dir,_score>(query, subject, mtx.score_buffer(), band, gap_open, gap_extend, j_best, i_best, max_score);
}

template<typename _val, typename _score, typename _traceback>
void floating_sw(const _val *query, local_match<_val> &segment, int band, _score xdrop, _score gap_open, _score gap_extend, const _traceback& = Score_only (), const _score& = int())
{
	segment += floating_sw_dir<_val,Right,_score,_traceback>(query, segment.subject_, band, xdrop, gap_open, gap_extend);
	const local_match<_val> left (floating_sw_dir<_val,Left,_score,_traceback>(query, segment.subject_, band, xdrop, gap_open, gap_extend));
	if(left.query_len_ > 0) {
		if(segment.transcript_ == 0)
			segment.transcript_ = new Edit_transcript;
		segment -= left;
		segment.query_begin_--;
		segment.subject_begin_--;
		const _val q = *query, s = mask_critical(*segment.subject_);
		segment.score_ -= score_matrix::get().letter_score(q, s);
		if(q == s)
			segment.identities_--;
		else
			segment.mismatches_--;
		segment.len_--;
		segment.subject_len_--;
		segment.query_len_--;
	}
}

#endif /* FLOATING_SW_H_ */

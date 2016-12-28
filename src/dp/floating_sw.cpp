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

#include <limits>
#include "../basic/match.h"
#include "../dp/scalar_dp_matrix.h"
#include "../util/direction.h"
#include "scalar_traceback.h"
#include "scalar_dp_matrix.h"
#include "../basic/score_matrix.h"
#include "../align/align.h"

template<typename _score, typename _traceback> TLS_PTR typename Score_buffer<_score,_traceback>::Type* Scalar_dp_matrix<_score,_traceback>::score_ptr = 0;
template<typename _score, typename _traceback> TLS_PTR Double_buffer<_score>* Scalar_dp_matrix<_score,_traceback>::hgap_ptr = 0;

template struct Scalar_dp_matrix<int, Traceback>;

template<typename _dir, typename _score>
Hsp_data get_traceback(const Letter *query,
	const Letter *subject,
	const Growing_buffer<_score> &scores,
	int band,
	_score gap_open,
	_score gap_extend,
	int i,
	int j,
	_score score,
	const Traceback&)
{
	return traceback<_dir, _score>(query, subject, scores, band, gap_open, gap_extend, i, j, score);
}

template<typename _dir, typename _score>
Hsp_data get_traceback(const Letter *query,
	const Letter *subject,
	const Double_buffer<_score> &scores,
	int band,
	_score gap_open,
	_score gap_extend,
	int i,
	int j,
	_score score,
	const Score_only&)
{
	return Hsp_data((int)score);
}

template<typename _dir, typename _score, typename _traceback>
Hsp_data floating_sw_dir(const Letter *query, const Letter* subject, int band, _score xdrop, _score gap_open, _score gap_extend, uint64_t &cell_updates)
{
	using std::max;

	_score max_score = 0, column_max = 0;
	int j = 0, i_max = -1, j_best = -1, i_best = -1;
	Scalar_dp_matrix<_score, _traceback> mtx(band);
	const Letter *x = query, *y = subject;

	while (*y != '\xff' && max_score - column_max < xdrop) {
		typename Scalar_dp_matrix<_score, _traceback>::Column_iterator it = mtx.column(j, i_max);
		if (get_dir(x, it.row(), _dir()) == '\xff')
			break;
		_score vgap = Scalar_dp_matrix<_score, _traceback>::minus_inf;
		if (get_dir(x, i_max + 1, _dir()) == '\xff') {
			column_max = std::numeric_limits<_score>::min();
		}
		else {
			++i_max;
			column_max += score_matrix(*y, get_dir(x, i_max, _dir()));
		}

		for (; it.valid() && get_dir(x, it.row(), _dir()) != '\xff'; ++it) {
			const _score match_score = (_score)score_matrix(*y, get_dir(x, it.row(), _dir()));
			const _score s = max(max(it.diag() + match_score, vgap), it.hgap_in());
			if (s > column_max) {
				column_max = s;
				i_max = it.row();
			}
			const _score open = s - gap_open;
			vgap = max(vgap - gap_extend, open);
			it.hgap_out() = max(it.hgap_in() - gap_extend, open);
			it.score() = s;
			++cell_updates;
		}

		if (column_max > max_score) {
			max_score = column_max;
			j_best = j;
			i_best = i_max;
		}
		y = inc_dir(y, _dir());
		++j;
	}
	
	return get_traceback<_dir, _score>(query, subject, mtx.score_buffer(), band, gap_open, gap_extend, j_best, i_best, max_score, _traceback());
}

template<typename _score, typename _traceback>
void floating_sw(const Letter *query, const Letter *subject, Hsp_data &segment, int band, _score xdrop, _score gap_open, _score gap_extend, uint64_t &cell_updates, unsigned query_anchor, unsigned subject_anchor, const _traceback&, const _score&)
{
	segment.merge(floating_sw_dir<Right, _score, _traceback>(query + 1, subject + 1, band, xdrop, gap_open, gap_extend, cell_updates),
		floating_sw_dir<Left, _score, _traceback>(query, subject, band, xdrop, gap_open, gap_extend, cell_updates), query_anchor, subject_anchor);
}

template void floating_sw<int, Traceback>(const Letter *query, const Letter *subject, Hsp_data &segment, int band, int xdrop, int gap_open, int gap_extend, uint64_t &cell_updates, unsigned query_anchor, unsigned subject_anchor, const Traceback&, const int&);
template void floating_sw<int, Score_only>(const Letter *query, const Letter *subject, Hsp_data &segment, int band, int xdrop, int gap_open, int gap_extend, uint64_t &cell_updates, unsigned query_anchor, unsigned subject_anchor, const Score_only&, const int&);
template void floating_sw<float, Traceback>(const Letter *query, const Letter *subject, Hsp_data &segment, int band, float xdrop, float gap_open, float gap_extend, uint64_t &cell_updates, unsigned query_anchor, unsigned subject_anchor, const Traceback&, const float&);
template void floating_sw<float, Score_only>(const Letter *query, const Letter *subject, Hsp_data &segment, int band, float xdrop, float gap_open, float gap_extend, uint64_t &cell_updates, unsigned query_anchor, unsigned subject_anchor, const Score_only&, const float&);
template void floating_sw<double, Traceback>(const Letter *query, const Letter *subject, Hsp_data &segment, int band, double xdrop, double gap_open, double gap_extend, uint64_t &cell_updates, unsigned query_anchor, unsigned subject_anchor, const Traceback&, const double&);
template void floating_sw<double, Score_only>(const Letter *query, const Letter *subject, Hsp_data &segment, int band, double xdrop, double gap_open, double gap_extend, uint64_t &cell_updates, unsigned query_anchor, unsigned subject_anchor, const Score_only&, const double&);
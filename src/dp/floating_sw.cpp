/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#include <limits>
#include "../basic/match.h"
#include "../dp/scalar_dp_matrix.h"
#include "../util/direction.h"
#include "scalar_traceback.h"
#include "scalar_dp_matrix.h"
#include "../basic/score_matrix.h"
#include "../align/align.h"
#include "dp.h"

template<typename _score, typename _traceback> TLS_PTR typename Score_buffer<_score, _traceback>::Type* Scalar_dp_matrix<_score, _traceback>::score_ptr = 0;
template<typename _score, typename _traceback> TLS_PTR Double_buffer<_score>* Scalar_dp_matrix<_score, _traceback>::hgap_ptr = 0;

template struct Scalar_dp_matrix<int, Traceback>;

template<typename _dir, typename _score, typename _score_correction>
Hsp_data get_traceback(const Letter *query,
	const Letter *subject,
	const Growing_buffer<_score> &scores,
	int band,
	_score gap_open,
	_score gap_extend,
	int i,
	int j,
	int query_anchor,
	_score score,
	const Traceback&,
	const _score_correction &score_correction)
{
	return traceback<_dir, _score, _score_correction>(query, subject, scores, band, gap_open, gap_extend, i, j, query_anchor, score, score_correction);
}

template<typename _dir, typename _score, typename _score_correction>
Hsp_data get_traceback(const Letter *query,
	const Letter *subject,
	const Double_buffer<_score> &scores,
	int band,
	_score gap_open,
	_score gap_extend,
	int i,
	int j,
	int query_anchor,
	_score score,
	const Score_only&,
	const _score_correction &score_correction)
{
	return Hsp_data((int)score);
}

template<typename _dir, typename _score, typename _traceback, typename _score_correction>
Hsp_data floating_sw_dir(const Letter *query, const Letter* subject, int band, _score xdrop, _score gap_open, _score gap_extend, uint64_t &cell_updates, const _score_correction &score_correction, int query_anchor, int min_j)
{
	using std::max;

	_score max_score = 0, column_max = 0;
	int j = 0, i_max = -1, j_best = -1, i_best = -1;
	Scalar_dp_matrix<_score, _traceback> mtx(band);
	const Letter *x = query, *y = subject;
	while (*y != '\xff' && (max_score - column_max < xdrop || j < min_j)) {
		typename Scalar_dp_matrix<_score, _traceback>::Column_iterator it = mtx.column(j, i_max);
		if (get_dir(x, it.row(), _dir()) == '\xff')
			break;
		_score vgap = Scalar_dp_matrix<_score, _traceback>::minus_inf;
		if (get_dir(x, i_max + 1, _dir()) == '\xff') {
			column_max = std::numeric_limits<_score>::min();
		}
		else {
			++i_max;
			_score match_score = (_score)score_matrix(*y, get_dir(x, i_max, _dir()));
			score_correction(match_score, i_max, query_anchor, _dir::mult);
			column_max += match_score;
		}

		for (; it.valid() && get_dir(x, it.row(), _dir()) != '\xff'; ++it) {
			_score match_score = (_score)score_matrix(*y, get_dir(x, it.row(), _dir()));
			score_correction(match_score, it.row(), query_anchor, _dir::mult);
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
	
	return get_traceback<_dir, _score, _score_correction>(query, subject, mtx.score_buffer(), band, gap_open, gap_extend, j_best, i_best, query_anchor, max_score, _traceback(), score_correction);
}

template<typename _score, typename _traceback, typename _score_correction>
void floating_sw(const Letter *query, const Letter *subject, Hsp_data &segment, int band, _score xdrop, _score gap_open, _score gap_extend, uint64_t &cell_updates, unsigned query_anchor, unsigned subject_anchor, int min_j, const _score_correction &score_correction, const _traceback&, const _score&)
{
	segment.merge(floating_sw_dir<Right, _score, _traceback, _score_correction>(query + 1, subject + 1, band, xdrop, gap_open, gap_extend, cell_updates, score_correction, query_anchor + 1, min_j-subject_anchor),
		floating_sw_dir<Left, _score, _traceback, _score_correction>(query, subject, band, xdrop, gap_open, gap_extend, cell_updates, score_correction, query_anchor,0), query_anchor, subject_anchor);
}

template void floating_sw<int, Traceback, No_score_correction>(const Letter *query, const Letter *subject, Hsp_data &segment, int band, int xdrop, int gap_open, int gap_extend, uint64_t &cell_updates, unsigned query_anchor, unsigned subject_anchor, int min_j, const No_score_correction&, const Traceback&, const int&);
template void floating_sw<int, Score_only, No_score_correction>(const Letter *query, const Letter *subject, Hsp_data &segment, int band, int xdrop, int gap_open, int gap_extend, uint64_t &cell_updates, unsigned query_anchor, unsigned subject_anchor, int min_j, const No_score_correction&, const Score_only&, const int&);
template void floating_sw<float, Traceback, Bias_correction>(const Letter *query, const Letter *subject, Hsp_data &segment, int band, float xdrop, float gap_open, float gap_extend, uint64_t &cell_updates, unsigned query_anchor, unsigned subject_anchor, int min_j, const Bias_correction&, const Traceback&, const float&);
template void floating_sw<float, Score_only, Bias_correction>(const Letter *query, const Letter *subject, Hsp_data &segment, int band, float xdrop, float gap_open, float gap_extend, uint64_t &cell_updates, unsigned query_anchor, unsigned subject_anchor, int min_j, const Bias_correction&, const Score_only&, const float&);
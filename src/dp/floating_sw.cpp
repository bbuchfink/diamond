#include <limits>
#include "../basic/match.h"
#include "../dp/scalar_dp_matrix.h"
#include "../util/direction.h"
#include "scalar_traceback.h"
#include "scalar_dp_matrix.h"
#include "../basic/score_matrix.h"

template<typename _score, typename _traceback> TLS_PTR typename Score_buffer<_score,_traceback>::Type* Scalar_dp_matrix<_score,_traceback>::score_ptr = 0;
template<typename _score, typename _traceback> TLS_PTR Double_buffer<_score>* Scalar_dp_matrix<_score,_traceback>::hgap_ptr = 0;

template struct Scalar_dp_matrix<int, Traceback>;

template<typename _dir, typename _score>
local_match get_traceback(const Letter *query,
	const Letter *subject,
	const Growing_buffer<_score> &scores,
	int band,
	int gap_open,
	int gap_extend,
	int i,
	int j,
	int score,
	vector<char> &transcript_buf,
	const Traceback&)
{
	return traceback<_dir, _score>(query, subject, scores, band, gap_open, gap_extend, i, j, score, transcript_buf);
}

template<typename _dir, typename _score>
local_match get_traceback(const Letter *query,
	const Letter *subject,
	const Double_buffer<_score> &scores,
	int band,
	int gap_open,
	int gap_extend,
	int i,
	int j,
	int score,
	vector<char> &transcript_buf,
	const Score_only&)
{
	return local_match(score);
}

template<typename _dir, typename _score, typename _traceback>
local_match floating_sw_dir(const Letter *query, const Letter* subject, int band, _score xdrop, _score gap_open, _score gap_extend, vector<char> &transcript_buf, unsigned long &cell_updates)
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
		_score vgap = Scalar_dp_matrix<_score, _traceback>::NEG_MIN;
		if (get_dir(x, i_max + 1, _dir()) == '\xff') {
			column_max = std::numeric_limits<_score>::min();
		}
		else {
			++i_max;
			column_max += score_matrix(mask_critical(*y), get_dir(x, i_max, _dir()));
		}

		for (; it.valid() && get_dir(x, it.row(), _dir()) != '\xff'; ++it) {
			const _score match_score = score_matrix(mask_critical(*y), get_dir(x, it.row(), _dir()));
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
	
	return get_traceback<_dir, _score>(query, subject, mtx.score_buffer(), band, gap_open, gap_extend, j_best, i_best, max_score, transcript_buf, _traceback());
}

template<typename _score, typename _traceback>
void floating_sw(const Letter *query, local_match &segment, int band, _score xdrop, _score gap_open, _score gap_extend, vector<char> &transcript_buf, unsigned long &cell_updates, const _traceback&, const _score&)
{
	segment += floating_sw_dir<Right, _score, _traceback>(query, segment.subject_, band, xdrop, gap_open, gap_extend, transcript_buf, cell_updates);
	const local_match left(floating_sw_dir<Left, _score, _traceback>(query, segment.subject_, band, xdrop, gap_open, gap_extend, transcript_buf, cell_updates));
	if (left.query_len_ > 0) {
		segment -= left;
		segment.query_begin_--;
		segment.subject_begin_--;
		const Letter q = *query, s = mask_critical(*segment.subject_);
		segment.score_ -= score_matrix(q, s);
		if (q == s)
			segment.identities_--;
		else
			segment.mismatches_--;
		segment.len_--;
		segment.subject_len_--;
		segment.query_len_--;
	}
}

template void floating_sw<int, Traceback>(const Letter *query, local_match &segment, int band, int xdrop, int gap_open, int gap_extend, vector<char> &transcript_buf, unsigned long &cell_updates, const Traceback&, const int&);
template void floating_sw<int, Score_only>(const Letter *query, local_match &segment, int band, int xdrop, int gap_open, int gap_extend, vector<char> &transcript_buf, unsigned long &cell_updates, const Score_only&, const int&);
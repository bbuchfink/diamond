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

#include "dp.h"
#include "../basic/score_matrix.h"

/*int xdrop_ungapped(const Letter *query, const Letter *subject, unsigned seed_len, unsigned &delta, unsigned &len)
{
	int score(0), st(0);
	unsigned n(0);
	delta = 0;

	const Letter *q(query - 1), *s(subject - 1);
	const unsigned window_left = std::max(config.window, (unsigned)config.seed_anchor) - config.seed_anchor;
	while (score - st < config.raw_ungapped_xdrop
		&& delta < window_left
		&& *q != sequence::DELIMITER
		&& *s != sequence::DELIMITER)
	{
#ifdef SEQ_MASK
		st += score_matrix(letter_mask(*q), letter_mask(*s));
#else
		st += score_matrix(*q, *s);
#endif
		score = std::max(score, st);
		--q;
		--s;
		++delta;
	}

	q = query + seed_len;
	s = subject + seed_len;
	st = score;
	assert(seed_len >= config.seed_anchor);
	const unsigned window_right = std::max(config.window, seed_len - config.seed_anchor) - (seed_len - config.seed_anchor);
	while (score - st < config.raw_ungapped_xdrop
		&& n < window_right
		&& *q != sequence::DELIMITER
		&& *s != sequence::DELIMITER)
	{
#ifdef SEQ_MASK
		st += score_matrix(letter_mask(*q), letter_mask(*s));
#else
		st += score_matrix(*q, *s);
#endif
		score = std::max(score, st);
		++q;
		++s;
		++n;
	}

	for (unsigned i = 0; i<seed_len; ++i)
#ifdef SEQ_MASK
		score += score_matrix(letter_mask(query[i]), letter_mask(subject[i]));
#else
		score += score_matrix(query[i], subject[i]);
#endif

	len = delta + n + seed_len;
	return score;
}*/

int xdrop_ungapped(const Letter *query, const Letter *subject, unsigned &delta, unsigned &len)
{
	int score(0), st(0), n=1;
	delta = 0;
	
	const Letter *q(query - 1), *s(subject - 1);
	while (score - st < config.raw_ungapped_xdrop
		&& *q != sequence::DELIMITER
		&& *s != sequence::DELIMITER)
	{
#ifdef SEQ_MASK
		st += score_matrix(letter_mask(*q), letter_mask(*s));
#else
		st += score_matrix(*q, *s);
#endif
		if (st > score) {
			score = st;
			delta = n;
		}
		--q;
		--s;
		++n;
	}

	q = query;
	s = subject;
	st = score;
	n = 1;
	len = 0;
	while (score - st < config.raw_ungapped_xdrop
		&& *q != sequence::DELIMITER
		&& *s != sequence::DELIMITER)
	{
#ifdef SEQ_MASK
		st += score_matrix(letter_mask(*q), letter_mask(*s));
#else
		st += score_matrix(*q, *s);
#endif
		if (st > score) {
			score = st;
			len = n;
		}
		++q;
		++s;
		++n;
	}
	len += delta;
	return score;
}

Diagonal_segment xdrop_ungapped(const sequence &query, const Bias_correction &query_bc, const sequence &subject, int qa, int sa)
{
	const float xdrop = (float)config.raw_ungapped_xdrop;
	float score = 0, st = 0;
	int n = 1, delta = 0, len = 0;

	int q = qa - 1, s = sa - 1;
	Letter ql, sl;
	while (score - st < xdrop
		&& (ql = query[q]) != sequence::DELIMITER
		&& (sl = subject[s]) != sequence::DELIMITER)
	{
		st += score_matrix(ql, sl) + query_bc[q];
		if (st > score) {
			score = st;
			delta = n;
		}
		--q;
		--s;
		++n;
	}

	q = qa;
	s = sa;
	st = score;
	n = 1;
	while (score - st < xdrop
		&& (ql = query[q]) != sequence::DELIMITER
		&& (sl = subject[s]) != sequence::DELIMITER)
	{
		st += score_matrix(ql, sl) + query_bc[q];
		if (st > score) {
			score = st;
			len = n;
		}
		++q;
		++s;
		++n;
	}
	return Diagonal_segment(qa - delta, sa - delta, len + delta, (int)score);
}

Diagonal_segment xdrop_ungapped(const sequence &query, const sequence &subject, int qa, int sa)
{
	const int xdrop = config.raw_ungapped_xdrop;
	int score = 0, st = 0;
	int n = 1, delta = 0, len = 0;

	int q = qa - 1, s = sa - 1;
	Letter ql, sl;
	while (score - st < xdrop
		&& (ql = query[q]) != sequence::DELIMITER
		&& (sl = subject[s]) != sequence::DELIMITER)
	{
		st += score_matrix(ql, sl);
		if (st > score) {
			score = st;
			delta = n;
		}
		--q;
		--s;
		++n;
	}

	q = qa;
	s = sa;
	st = score;
	n = 1;
	while (score - st < xdrop
		&& (ql = query[q]) != sequence::DELIMITER
		&& (sl = subject[s]) != sequence::DELIMITER)
	{
		st += score_matrix(ql, sl);
		if (st > score) {
			score = st;
			len = n;
		}
		++q;
		++s;
		++n;
	}
	return Diagonal_segment(qa - delta, sa - delta, len + delta, score);
}

int xdrop_ungapped_right(const Letter *query, const Letter *subject, int &len)
{
	int score(0), st(0), n = 1;
	len = 0;
	
	const Letter *q = query;
	const Letter *s = subject;
	
	while (score - st < config.raw_ungapped_xdrop
		&& *q != sequence::DELIMITER
		&& *s != sequence::DELIMITER)
	{
#ifdef SEQ_MASK
		st += score_matrix(letter_mask(*q), letter_mask(*s));
#else
		st += score_matrix(*q, *s);
#endif
		if (st > score) {
			score = st;
			len = n;
		}
		++q;
		++s;
		++n;
	}
	return score;
}

int ungapped_window(const Letter* query, const Letter* subject, int window) {
	int score = 0, st = 0, n = 0;
	const Letter* q = query, * s = subject;
	while (n < window)
	{
#ifdef SEQ_MASK
		st += score_matrix(letter_mask(*q), *s);
#else
		st += score_matrix(*q, *s);
#endif
		st = std::max(st, 0);
		score = std::max(score, st);
		++q;
		++s;
		++n;
	}
	return score;
}

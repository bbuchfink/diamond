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

#include <cstddef>
#include "dp.h"
#include "../stats/score_matrix.h"
#include "../stats/hauser_correction.h"
#include "ungapped.h"
#include "../util/sequence/sequence.h"

using std::max;
using std::pair;
using std::min;
using std::nullptr_t;
using std::vector;

int xdrop_ungapped(const Letter *query, const Letter *subject, unsigned &delta, unsigned &len)
{
	int score(0), st(0), n=1;
	delta = 0;
	
	const Letter *q(query - 1), *s(subject - 1);
	while (score - st < config.raw_ungapped_xdrop
		&& *q != Sequence::DELIMITER
		&& *s != Sequence::DELIMITER)
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
		&& *q != Sequence::DELIMITER
		&& *s != Sequence::DELIMITER)
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

DiagonalSegment xdrop_ungapped(const Sequence &query, const Bias_correction &query_bc, const Sequence &subject, int qa, int sa)
{
	const float xdrop = (float)config.raw_ungapped_xdrop;
	float score = 0, st = 0;
	int n = 1, delta = 0, len = 0;

	int q = qa - 1, s = sa - 1;
	Letter ql, sl;
	while (score - st < xdrop
		&& (ql = query[q]) != Sequence::DELIMITER
		&& (sl = subject[s]) != Sequence::DELIMITER)
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
		&& (ql = query[q]) != Sequence::DELIMITER
		&& (sl = subject[s]) != Sequence::DELIMITER)
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
	return DiagonalSegment(qa - delta, sa - delta, len + delta, (int)score);
}

static void add_cbs(Score&, Loc, nullptr_t) {}

static void add_cbs(Score& score, Loc loc, const int8_t* cbs) {
	score += (Score)cbs[loc];
}

struct CountIdentities {
	CountIdentities() :
		n(0)
	{}
	void operator()(Letter a, Letter b) {
		const bool eq = a == b;
		n += (Loc)eq;
	}
	void update(Loc& ident) {
		ident += n;
		n = 0;
	}
	Loc n;
};

struct ScoreOnly {
	void operator()(Letter a, Letter b) const {}
	void update(Loc) const {}
};

template<typename Id, typename Cbs>
DiagonalSegment xdrop_ungapped(const Sequence& query, Cbs query_cbs, const Sequence& subject, int qa, int sa, const Id&)
{
	Id id;
	const int xdrop = config.raw_ungapped_xdrop;
	int score = 0, st = 0;
	int n = 1, delta = 0, len = 0, ident = 0;

	int q = qa - 1, s = sa - 1;
	Letter ql, sl;
	while (score - st < xdrop
		&& (ql = query[q]) != Sequence::DELIMITER
		&& (sl = subject[s]) != Sequence::DELIMITER)
	{
		st += score_matrix(ql, sl);
		add_cbs(st, q, query_cbs);
		id(ql, sl);
		if (st > score) {
			score = st;
			delta = n;
			id.update(ident);
		}
		--q;
		--s;
		++n;
	}

	q = qa;
	s = sa;
	st = score;
	n = 1;
	id = Id();
	while (score - st < xdrop
		&& (ql = query[q]) != Sequence::DELIMITER
		&& (sl = subject[s]) != Sequence::DELIMITER)
	{
		st += score_matrix(ql, sl);
		add_cbs(st, q, query_cbs);
		id(ql, sl);
		if (st > score) {
			score = st;
			len = n;
			id.update(ident);
		}
		++q;
		++s;
		++n;
	}
	return DiagonalSegment(qa - delta, sa - delta, len + delta, score, ident);
}

DiagonalSegment xdrop_ungapped(const Sequence& query, const int8_t* query_cbs, const Sequence& subject, int qa, int sa, bool count_identities) {
	if (count_identities) {
		if (query_cbs == nullptr)
			return xdrop_ungapped(query, nullptr, subject, qa, sa, CountIdentities());
		else
			return xdrop_ungapped(query, query_cbs, subject, qa, sa, CountIdentities());
	}
	else {
		if (query_cbs == nullptr)
			return xdrop_ungapped(query, nullptr, subject, qa, sa, ScoreOnly());
		else
			return xdrop_ungapped(query, query_cbs, subject, qa, sa, ScoreOnly());
	}
}

int xdrop_ungapped_right(const Letter *query, const Letter *subject, int &len)
{
	int score(0), st(0), n = 1;
	len = 0;
	
	const Letter *q = query;
	const Letter *s = subject;
	
	while (score - st < config.raw_ungapped_xdrop
		&& *q != Sequence::DELIMITER
		&& *s != Sequence::DELIMITER)
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
		st += score_matrix(letter_mask(*q), letter_mask(*s));
		st = std::max(st, 0);
		score = std::max(score, st);
		++q;
		++s;
		++n;
	}
	return score;
}

Score self_score(const Sequence& seq)
{
	Score s = 0, sl = 0;
	if (Stats::CBS::hauser(config.comp_based_stats)) {
		Bias_correction cbs(seq);

		for (Loc i = 0; i < seq.length(); ++i) {
			const Letter l = seq[i];
			sl += score_matrix(l, l) + cbs.int8[i];
			sl = max(sl, 0);
			s = max(s, sl);
		}
	}
	else {
		for (Loc i = 0; i < seq.length(); ++i) {
			const Letter l = seq[i];
			sl += score_matrix(l, l);
			sl = max(sl, 0);
			s = max(s, sl);
		}
	}
	return s;
}

int score_range(Sequence query, Sequence subject, int i, int j, int j_end)
{
	int score = 0;
	while (j < j_end) {
		score += score_matrix(query[i], subject[j]);
		++i;
		++j;
	}
	return score;
}

template<typename Cbs>
DiagonalSegment score_range_s(Sequence query, Cbs query_cbs, Sequence subject, int i_begin, int j_begin, int j_end) {
	int score = 0, i = i_begin;
	for (Loc j = j_begin; j < j_end; ++j, ++i) {
		score += score_matrix(query[i], subject[j]);
		add_cbs(score, i, query_cbs);
	}
	return DiagonalSegment(i_begin, j_begin, j_end - j_begin, score);
}

template DiagonalSegment score_range_s<const int8_t*>(Sequence, const int8_t*, Sequence, int, int, int);
template DiagonalSegment score_range_s<nullptr_t>(Sequence, nullptr_t, Sequence, int, int, int);

template<typename T>
struct Increment {
	T operator()(T x) const {
		return ++x;
	}
};

template<typename T>
struct Decrement {
	T operator()(T x) const {
		return --x;
	}
};

template<typename Inc>
pair<Score, Loc> xdrop_anchored(const Letter* p1, const Letter* p2, Loc len, Inc inc) {
	if (len == 0)
		return { 0,0 };
	Score max_score = 0, score = 0;
	Loc max_n = 0, n = 0;
	do {
		score += score_matrix(letter_mask(*p1), letter_mask(*p2));
		++n;
		p1 = inc(p1);
		p2 = inc(p2);
		if (score > max_score) {
			max_score = score;
			max_n = n;
		}
	} while (n < len && max_score - score < config.raw_ungapped_xdrop);
	return { max_score, max_n };
}

DiagonalSegment xdrop_ungapped(const Sequence& query, const Sequence& subject, const DiagonalSegment& anchor) {
	auto left = xdrop_anchored(query.data() + anchor.i - 1, subject.data() + anchor.j - 1, std::min(anchor.i, anchor.j), Decrement<const Letter*>());
	auto right = xdrop_anchored(query.data() + anchor.query_end(),
		subject.data() + anchor.subject_end(),
		std::min(query.length() - anchor.query_end(), subject.length() - anchor.subject_end()),
		Increment<const Letter*>());
	return { anchor.i - left.second, anchor.j - left.second,
		anchor.len + left.second + right.second,
		left.first + right.first + score_range(query, subject, anchor.i, anchor.j, anchor.subject_end()) };
}

static Hsp trivial(Sequence query, Sequence target, Loc dq, Loc dt, const int8_t* query_cbs) {
	static const Loc WINDOW = 40, ID = 30;
	const Loc l = min(query.length() - dq, target.length() - dt);
	const uint64_t bits = ((uint64_t)1 << WINDOW) - 1;
	Loc n = 0;
	Score score = 0;
	uint64_t mask = 0;
	for (Loc i = 0; i < l; ++i) {
		uint64_t eq = query[i + dq] == target[i + dt];
		mask = ((mask << 1) | eq) & bits;
		++n;
		if (n >= WINDOW && popcount64(mask) < ID)
			return Hsp();
		score += score_matrix(query[i + dq], target[i + dt]);
		if (query_cbs)
			score += query_cbs[i + dq];
	}
	const double evalue = score_matrix.evalue(score, query.length(), target.length());
	if (evalue > config.max_evalue)
		return Hsp();
	if (l < WINDOW && (double)popcount64(mask) / l < ID / WINDOW)
		return Hsp();
	Hsp hsp;
	hsp.score = score;
	hsp.query_range = hsp.query_source_range = { dq, dq + l };
	hsp.subject_range = { dt, dt + l };
	hsp.evalue = evalue;
	hsp.bit_score = score_matrix.bitscore(score);
	return hsp;
}

Hsp trivial(Sequence query, Sequence target, const int8_t* query_cbs) {
	if (query.length() <= target.length()) {
		for (Loc i = 0; i <= target.length() - query.length(); ++i) {
			Hsp hsp = trivial(query, target, 0, i, query_cbs);
			if (hsp.score)
				return hsp;
		}
	}
	else {
		for (Loc i = 0; i <= query.length() - target.length(); ++i) {
			Hsp hsp = trivial(query, target, i, 0, query_cbs);
			if (hsp.score)
				return hsp;
		}
	}
	return Hsp();
}

Anchor make_clipped_anchor(const Anchor& anchor, Sequence query, const int8_t* query_cbs, Sequence target) {
	const Sequence q = query.subseq(anchor.query_begin(), anchor.query_end()),
		t = target.subseq(anchor.subject_begin(), anchor.subject_end());
	const vector<Score> s = Util::Seq::window_scores(q, t, config.anchor_window);
	const Score cutoff = (Score)round(config.anchor_score * config.anchor_window);
	const auto pred = [cutoff](Score s) { return s >= cutoff; };

	const auto max_window = max_element(s.begin(), s.end());
	const auto right = find_if_not(max_window + 1, s.cend(), pred);

	const auto max_window_r = s.crbegin() + (s.cend() - max_window - 1);
	const auto left_r = find_if_not(max_window_r + 1, s.crend(), pred);
	Loc d1 = Loc(right - s.cbegin()),
		d0 = std::max(Loc(s.crend() - left_r - config.anchor_window + 1), 0);
	while (d0 < q.length() && q[d0] != t[d0]) ++d0;
	while (d1 > 0 && q[d1 - 1] != t[d1 - 1]) --d1;
	if (d1 <= d0)
		return DiagonalSegment();
	const DiagonalSegment clipped_anchor = query_cbs == nullptr
		? score_range_s(query, nullptr, target, anchor.query_begin() + d0, anchor.subject_begin() + d0, anchor.subject_begin() + d1)
		: score_range_s(query, query_cbs, target, anchor.query_begin() + d0, anchor.subject_begin() + d0, anchor.subject_begin() + d1);
	const Score clipped_score = query_cbs == nullptr
		? score_range_s(query, nullptr, target, anchor.query_begin() + d1, anchor.subject_begin() + d1, anchor.subject_end()).score
		: score_range_s(query, query_cbs, target, anchor.query_begin() + d1, anchor.subject_begin() + d1, anchor.subject_end()).score;
	return Anchor(clipped_anchor, anchor.d_min_left, anchor.d_max_left, anchor.d_min_right, anchor.d_max_right, anchor.prefix_score - clipped_score);
}

Anchor make_null_anchor(const Anchor& anchor) {
	return Anchor(anchor.i + anchor.len / 2, anchor.j + anchor.len / 2, 0, 0);
}
/****
DIAMOND protein aligner
Copyright (C) 2016-2021 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
						
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
#include "../../basic/translated_position.h"
#include "../../basic/sequence.h"

struct DiagonalSegment
{
	DiagonalSegment() :
		len(0),
		score(0)
	{}
	DiagonalSegment(int query_pos, int subject_pos, int len, int score, Loc ident = 0) :
		i(query_pos),
		j(subject_pos),
		len(len),
		score(score),
		ident(ident)
	{}
	bool empty() const
	{
		return len == 0;
	}
	Interval query_range() const
	{
		return Interval(i, i + len);
	}
	Interval subject_range() const
	{
		return Interval(j, j + len);
	}
	int subject_begin() const {
		return j;
	}
	int subject_last() const
	{
		return j + len - 1;
	}
	int query_last() const
	{
		return i + len - 1;
	}
	int subject_end() const
	{
		return j + len;
	}
	int query_begin() const {
		return i;
	}
	int query_end() const
	{
		return i + len;
	}
	int diag() const
	{
		return i - j;
	}
	double id_percent() const {
		return (double)ident / len * 100.0;
	}
	double cov_percent(Loc seq_len) const {
		return (double)len / seq_len * 100.0;
	}
	Interval band_interval(const Loc band) const {
		return { diag() - band, diag() + band };
	}
	void set_query_end(Loc i) {
		len = i - this->i;
	}
	void set_target_end(Loc j) {
		len = j - this->j;
	}
	DiagonalSegment intersect(const DiagonalSegment&x) const
	{
		if (diag() != x.diag())
			return DiagonalSegment();
		else {
			const Interval q = ::intersect(query_range(), x.query_range());
			return DiagonalSegment(q.begin_, ::intersect(subject_range(), x.subject_range()).begin_, q.length(), 0);
		}
	}
	bool is_enveloped(const DiagonalSegment&x) const
	{
		return score <= x.score
			&& query_range().overlap_factor(x.query_range()) == 1
			&& subject_range().overlap_factor(x.subject_range()) == 1;
	}
	DiagonalSegment transpose() const
	{
		return DiagonalSegment(j, i, len, score);
	}
	int partial_score(int diff) const
	{
		return score*std::max(len - diff, 0) / len;
	}
	bool operator<=(const DiagonalSegment&rhs) const
	{
		return i + len <= rhs.i && j + len <= rhs.j;
	}
	bool operator==(const DiagonalSegment&rhs) const
	{
		return i == rhs.i && j == rhs.j && len == rhs.len;
	}
	static bool cmp_subject(const DiagonalSegment&x, const DiagonalSegment&y)
	{
		return x.j < y.j || (x.j == y.j && x.i < y.i);
	}
	static bool cmp_score(const DiagonalSegment& x, const DiagonalSegment& y)
	{
		return x.score > y.score;
	}
	static bool cmp_subject_end(const DiagonalSegment&x, const DiagonalSegment&y)
	{
		return x.subject_end() < y.subject_end();
	}
	static bool cmp_heuristic(const DiagonalSegment&x, const DiagonalSegment&y)
	{
		return (x.subject_end() < y.subject_end() && x.j < y.j)
			|| (x.j - y.j < y.subject_end() - x.subject_end());
	}
	static bool cmp_diag(const DiagonalSegment&x, const DiagonalSegment& y)
	{
		return x.diag() < y.diag() || (x.diag() == y.diag() && x.j < y.j);
	}
	static bool cmp_len(const DiagonalSegment& x, const DiagonalSegment& y) {
		return x.len > y.len;
	}
	friend int abs_shift(const DiagonalSegment&x, const DiagonalSegment& y)
	{
		return abs(x.diag() - y.diag());
	}
	friend std::ostream& operator<<(std::ostream &s, const DiagonalSegment& d)
	{
		s << "i=" << d.i << " j=" << d.j << " l=" << d.len << " score=" << d.score;
		return s;
	}
	int i, j, len, score, ident;
};

struct DiagonalSegmentT
{
	DiagonalSegmentT():
		len(0)
	{}

	DiagonalSegmentT(const TranslatedPosition &i, int j, int len, int score=0):
		i(i),
		j(j),
		len(len),
		score(score)
	{}

	DiagonalSegmentT(const DiagonalSegment& d, Frame frame):
		i(TranslatedPosition(d.i, frame)),
		j(d.j),
		len(d.len),
		score(d.score)
	{}

	int subject_last() const
	{
		return j + len - 1;
	}

	TranslatedPosition query_last() const
	{
		return i + len - 1;
	}

	int subject_end() const
	{
		return j + len;
	}

	TranslatedPosition query_end() const
	{
		return i + len;
	}

	int diag() const
	{
		return i - j;
	}

	friend std::ostream& operator<<(std::ostream &s, const DiagonalSegmentT &d)
	{
		s << "i=(" << d.i << ") j=" << d.j << " len=" << d.len << " score=" << d.score << std::endl;
		return s;
	}

	Interval query_absolute_range(int dna_len) const
	{
		return TranslatedPosition::absolute_interval(i, i + len, dna_len);
	}

	Interval query_in_strand_range() const
	{
		return Interval(i.in_strand(), (i + len).in_strand());
	}

	Interval subject_range() const
	{
		return Interval(j, j + len);
	}

	int partial_score(const DiagonalSegmentT &d) const
	{
		const double overlap = std::max(subject_range().overlap_factor(d.subject_range()), query_in_strand_range().overlap_factor(d.query_in_strand_range()));
		return int((1.0 - overlap)*score);
	}

	void cut_out(const DiagonalSegmentT &d)
	{
		const int ll = std::min(d.i.translated - i.translated, d.j - j),
			lr = std::min(query_end().translated - d.query_end().translated, subject_end() - d.subject_end());
		int len2;
		if (ll > 0 && ll >= lr) {
			len2 = std::min(len, ll);
		}
		else if (lr > 0 && lr >= ll) {
			len2 = std::min(len, lr);
			i = query_end() - len2;
			j = subject_end() - len2;
		}
		else {
			len2 = 0;
		}
		score = int((double)len2 / len * score);
		len = len2;
	}

	TranslatedPosition i;
	int j, len, score;
};

inline int diag_count(Loc query_len, Loc target_len) {
	return query_len + target_len - 1;
}

inline int diag_idx(Loc i, Loc j, Loc target_len) {
	return i - j + target_len - 1;
}
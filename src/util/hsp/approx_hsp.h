/****
DIAMOND protein aligner
Copyright (C) 2016-2022 Max Planck Society for the Advancement of Science e.V.
						Benjamin Buchfink

Code developed by Benjamin Buchfink <buchfink@gmail.com>

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
#include <limits.h>
#include <algorithm>
#include <float.h>
#include "../util/geo/interval.h"
#include "../util/geo/diagonal_segment.h"

struct Anchor : public DiagonalSegment {
	Anchor():
		DiagonalSegment(),
		prefix_score(0),
		d_min_left(std::numeric_limits<Loc>::max()), d_max_left(std::numeric_limits<Loc>::min()), d_min_right(std::numeric_limits<Loc>::max()), d_max_right(std::numeric_limits<Loc>::min())
	{}
	Anchor(const DiagonalSegment& d) :
		DiagonalSegment(d),
		prefix_score(0),
		d_min_left(std::numeric_limits<Loc>::max()), d_max_left(std::numeric_limits<Loc>::min()), d_min_right(std::numeric_limits<Loc>::max()), d_max_right(std::numeric_limits<Loc>::min())
	{
	}
	Anchor(const DiagonalSegment& d, Loc d_min_left, Loc d_max_left, Loc d_min_right, Loc d_max_right, Score prefix_score) :
		DiagonalSegment(d),
		prefix_score(prefix_score),
		d_min_left(d_min_left), d_max_left(d_max_left), d_min_right(d_min_right), d_max_right(d_max_right)
	{
	}
	Anchor(int query_pos, int subject_pos, int len, int score, Loc ident = 0) :
		DiagonalSegment(query_pos, subject_pos, len, score, ident),
		prefix_score(0),
		d_min_left(std::numeric_limits<Loc>::max()), d_max_left(std::numeric_limits<Loc>::min()), d_min_right(std::numeric_limits<Loc>::max()), d_max_right(std::numeric_limits<Loc>::min())
	{
	}
	Anchor& operator=(const DiagonalSegment& d) {
		*(DiagonalSegment*)this = d;
		return *this;
	}
	Score prefix_score;
	Loc d_min_left, d_max_left, d_min_right, d_max_right;
};

struct ApproxHsp
{
	ApproxHsp(unsigned frame, Score score = 0) :
		d_min(std::numeric_limits<int>::max()),
		d_max(std::numeric_limits<int>::min()),
		score(score),
		frame((int)frame)
	{}
	ApproxHsp(const Interval &query_source_range) :
		query_source_range(query_source_range)
	{}
	ApproxHsp(int d_min, int d_max, int score, int frame, const Interval& query_range, const Interval& subject_range, const Anchor& max_diag, double evalue = DBL_MAX):
		d_min(d_min),
		d_max(d_max),
		score(score),
		frame(frame),
		query_range(query_range),
		subject_range(subject_range),
		evalue(evalue),
		max_diag(max_diag)
	{}
	int partial_score(const DiagonalSegment &d) const
	{
		const double overlap = std::max(d.subject_range().overlap_factor(subject_range), d.query_range().overlap_factor(query_range));
		return int((1 - overlap)*d.score);
	}
	int partial_score(const ApproxHsp &x) const
	{
		const double overlap = std::max(x.subject_range.overlap_factor(subject_range), x.query_range.overlap_factor(query_range));
		return int((1 - overlap)*x.score);
	}
	bool disjoint(const DiagonalSegment &d) const
	{
		return intersect(query_range, d.query_range()).length() == 0 && intersect(subject_range, d.subject_range()).length() == 0;
	}
	bool disjoint(const ApproxHsp &x) const
	{
		return intersect(query_range, x.query_range).length() == 0 && intersect(subject_range, x.subject_range).length() == 0;
	}
	bool rel_disjoint(const DiagonalSegment &d) const
	{
		return intersect(query_range, d.query_range()).length() == 0 || intersect(subject_range, d.subject_range()).length() == 0;
	}
	bool rel_disjoint(const ApproxHsp &x) const
	{
		return intersect(query_range, x.query_range).length() == 0 || intersect(subject_range, x.subject_range).length() == 0;
	}
	bool collinear(const ApproxHsp &x) const
	{
		const int di = x.query_range.begin_ - query_range.begin_, dj = x.subject_range.begin_ - subject_range.begin_;
		return (di >= 0 && dj >= 0) || (di <= 0 && dj <= 0);
	}
	bool collinear(const DiagonalSegment &d) const
	{
		const int di = d.i - query_range.begin_, dj = d.j - subject_range.begin_;
		return (di >= 0 && dj >= 0) || (di <= 0 && dj <= 0);
	}
	static bool cmp_diag(const ApproxHsp& x, const ApproxHsp& y)
	{
		return x.frame < y.frame || (x.frame == y.frame && x.d_min < y.d_min);
	}
	friend std::ostream& operator<<(std::ostream& s, const ApproxHsp& h) {
		s << "Score=" << h.score << " query_range=" << h.query_range << " target_range=" << h.subject_range << std::endl;
		return s;
	}
	struct Frame
	{
		unsigned operator()(const ApproxHsp &x) const
		{
			return x.frame;
		}
	};
	int d_min, d_max, score, frame;
	Interval query_source_range, query_range, subject_range;
	double evalue;
	Anchor max_diag;
};

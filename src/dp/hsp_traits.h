#ifndef HSP_TRAITS_H_
#define HSP_TRAITS_H_

#include <limits.h>
#include <algorithm>
#include "../util/interval.h"
#include "../basic/diagonal_segment.h"
#include "../basic/match.h"

struct Hsp_traits
{
	Hsp_traits(unsigned frame) :
		d_min(std::numeric_limits<int>::max()),
		d_max(std::numeric_limits<int>::min()),
		score(0),
		frame((int)frame)
	{}
	Hsp_traits(const interval &query_source_range) :
		query_source_range(query_source_range)
	{}
	Hsp_traits(const Hsp &hsp)
	{}
	Hsp_traits(int d_min, int d_max, int score, int frame, const interval& query_range, const interval& subject_range):
		d_min(d_min),
		d_max(d_max),
		score(score),
		frame(frame),
		query_range(query_range),
		subject_range(subject_range)
	{}
	int partial_score(const Diagonal_segment &d) const
	{
		const double overlap = std::max(d.subject_range().overlap_factor(subject_range), d.query_range().overlap_factor(query_range));
		return int((1 - overlap)*d.score);
	}
	int partial_score(const Hsp_traits &x) const
	{
		const double overlap = std::max(x.subject_range.overlap_factor(subject_range), x.query_range.overlap_factor(query_range));
		return int((1 - overlap)*x.score);
	}
	bool disjoint(const Diagonal_segment &d) const
	{
		return intersect(query_range, d.query_range()).length() == 0 && intersect(subject_range, d.subject_range()).length() == 0;
	}
	bool disjoint(const Hsp_traits &x) const
	{
		return intersect(query_range, x.query_range).length() == 0 && intersect(subject_range, x.subject_range).length() == 0;
	}
	bool rel_disjoint(const Diagonal_segment &d) const
	{
		return intersect(query_range, d.query_range()).length() == 0 || intersect(subject_range, d.subject_range()).length() == 0;
	}
	bool rel_disjoint(const Hsp_traits &x) const
	{
		return intersect(query_range, x.query_range).length() == 0 || intersect(subject_range, x.subject_range).length() == 0;
	}
	bool collinear(const Hsp_traits &x) const
	{
		const int di = x.query_range.begin_ - query_range.begin_, dj = x.subject_range.begin_ - subject_range.begin_;
		return (di >= 0 && dj >= 0) || (di <= 0 && dj <= 0);
	}
	bool collinear(const Diagonal_segment &d) const
	{
		const int di = d.i - query_range.begin_, dj = d.j - subject_range.begin_;
		return (di >= 0 && dj >= 0) || (di <= 0 && dj <= 0);
	}
	static bool cmp_diag(const Hsp_traits& x, const Hsp_traits& y)
	{
		return x.frame < y.frame || (x.frame == y.frame && x.d_min < y.d_min);
	}
	struct Frame
	{
		unsigned operator()(const Hsp_traits &x) const
		{
			return x.frame;
		}
	};
	int d_min, d_max, score, frame;
	interval query_source_range, query_range, subject_range;
};

#endif
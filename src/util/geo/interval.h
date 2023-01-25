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

#pragma once
#include <ostream>
#include <stdexcept>
#include <algorithm>
#include <string>
#include "../../basic/value.h"

struct Interval
{
	Interval() :
		begin_(0),
		end_(0)
	{ }
	Interval(int begin, int end) :
		begin_(begin),
		end_(end)
	{ }
	int length() const
	{
		return end_ > begin_ ? end_ - begin_ : 0;
	}
	unsigned overlap(const Interval &rhs) const
	{
		return intersect(*this, rhs).length();
	}
	double overlap_factor(const Interval &rhs) const
	{
		return (double)overlap(rhs) / (double)length();
	}
	bool includes(int p) const
	{
		return p >= begin_ && p < end_;
	}
	bool contains(const Interval &i) const
	{
		return begin_ <= i.begin_ && end_ >= i.end_;
	}
	friend std::ostream& operator<<(std::ostream &os, const Interval &x)
	{
		os << "[" << x.begin_ << ";" << x.end_ << "]"; return os;
	}
	bool operator<(const Interval &rhs) const
	{
		return begin_ < rhs.begin_;
	}
	void merge(const Interval &k)
	{
		begin_ = std::min(begin_, k.begin_);
		end_ = std::max(end_, k.end_);
	}
	void check(int len) const {
		if (begin_ < 0 || end_ < 0 || end_ < begin_ || begin_ >= len || end_ > len)
			throw std::out_of_range("");
	}
	friend Interval intersect(const Interval &lhs, const Interval &rhs);
	int begin_, end_;
};

inline Interval intersect(const Interval &lhs, const Interval &rhs)
{
	return Interval(std::max(lhs.begin_, rhs.begin_), std::min(lhs.end_, rhs.end_));
}

template<typename It, typename Out>
void make_disjoint(const It begin, const It end, Out out) {
	if (begin == end)
		return;
	std::sort(begin, end);
	Loc a = begin->begin_, b = begin->end_;
	for (It i = std::next(begin); i != end; ++i) {
		if (i->begin_ <= b)
			b = std::max(b, i->end_);
		else {
			*out++ = Interval{ a,b };
			a = i->begin_;
			b = i->end_;
		}
	}
	*out = Interval{ a,b };
}
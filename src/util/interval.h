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

#ifndef INTERVAL_H_
#define INTERVAL_H_

#include <iostream>
#include <algorithm>

struct interval
{
	interval() :
		begin_(0),
		end_(0)
	{ }
	interval(int begin, int end) :
		begin_(begin),
		end_(end)
	{ }
	int length() const
	{
		return end_ > begin_ ? end_ - begin_ : 0;
	}
	unsigned overlap(const interval &rhs) const
	{
		return intersect(*this, rhs).length();
	}
	double overlap_factor(const interval &rhs) const
	{
		return (double)overlap(rhs) / (double)length();
	}
	bool includes(int p) const
	{
		return p >= begin_ && p < end_;
	}
	bool contains(const interval &i) const
	{
		return begin_ <= i.begin_ && end_ >= i.end_;
	}
	friend std::ostream& operator<<(std::ostream &os, const interval &x)
	{
		os << "[" << x.begin_ << ";" << x.end_ << "]"; return os;
	}
	bool operator<(const interval &rhs) const
	{
		return begin_ < rhs.begin_;
	}
	void merge(const interval &k)
	{
		begin_ = std::min(begin_, k.begin_);
		end_ = std::max(end_, k.end_);
	}
	friend interval intersect(const interval &lhs, const interval &rhs);
	int begin_, end_;
};

inline interval intersect(const interval &lhs, const interval &rhs)
{
	return interval(std::max(lhs.begin_, rhs.begin_), std::min(lhs.end_, rhs.end_));
}

#endif
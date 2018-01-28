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

#ifndef INTERVAL_PARTITION_H_
#define INTERVAL_PARTITION_H_

#include <assert.h>
#include <map>
#include <limits.h>
#include "interval.h"

struct IntervalNode
{
	IntervalNode():
		count(0),
		min_score(INT_MAX),
		max_score(0)
	{}
	IntervalNode(int count, int min_score, int max_score):
		count(count),
		min_score(min_score),
		max_score(max_score)
	{
	}
	IntervalNode add(int score, int cap) const
	{
		return IntervalNode(count + 1, count < cap ? std::min(min_score, score) : min_score, std::max(max_score, score));
	}
	int count, min_score, max_score;
};

struct IntervalPartition : protected std::map<int, IntervalNode>
{

	struct MaxScore {};
	struct MinScore {};

	struct Iterator
	{
		Iterator(const_iterator i, const_iterator j, const IntervalPartition &parent) :
			i_(i),
			j_(j),
			parent_(parent)
		{
		}
		bool good() const
		{
			return i_ != parent_.end();
		}
		Iterator& operator++()
		{
			i_ = j_;
			++j_;
			return *this;
		}
		std::pair<interval, IntervalNode> operator*() const
		{
			return std::make_pair(interval(i_->first, j_ == parent_.end() ? INT_MAX : j_->first), i_->second);
		}
	private:
		const_iterator i_, j_;
		const IntervalPartition &parent_;
	};

	IntervalPartition(int cap):
		cap(cap)
	{
		(*this)[0] = IntervalNode();
	}

	void insert(interval k, int score)
	{
		iterator i = lower_bound(k.begin_);
		if (i == end())
			i = std::map<int, IntervalNode>::insert(std::make_pair(k.begin_, IntervalNode())).first;
		else if (i->first != k.begin_) {
			//assert(i != std::map<int, IntervalNode>::begin());
			i--;
			i = std::map<int, IntervalNode>::insert(std::make_pair(k.begin_, i->second)).first;
		}
		IntervalNode last;
		while (i != end() && i->first < k.end_) {
			last = i->second;
			i->second = i->second.add(score, cap);
			++i;
		}
		if (i == end() || i->first != k.end_)
			(*this)[k.end_] = last;
	}

	int covered(interval k) const
	{
		Iterator i = begin(k.begin_);
		std::pair<interval, IntervalNode> l;
		int c = 0;
		while (i.good() && (l = *i).first.begin_ < k.end_) {
			if (l.second.count >= cap)
				c += k.overlap(l.first);
			++i;
		}
		return c;
	}

	int covered(interval k, int max_score, const MaxScore&) const
	{
		Iterator i = begin(k.begin_);
		std::pair<interval, IntervalNode> l;
		int c = 0;
		while (i.good() && (l = *i).first.begin_ < k.end_) {
			if (l.second.max_score >= max_score)
				c += k.overlap(l.first);
			++i;
		}
		return c;
	}

	int covered(interval k, int min_score, const MinScore&) const
	{
		Iterator i = begin(k.begin_);
		std::pair<interval, IntervalNode> l;
		int c = 0;
		while (i.good() && (l = *i).first.begin_ < k.end_) {
			if (l.second.count >= cap && l.second.min_score >= min_score)
				c += k.overlap(l.first);
			++i;
		}
		return c;
	}

	int min_score(interval k) const
	{
		Iterator i = begin(k.begin_);
		std::pair<interval, IntervalNode> l;
		int s = INT_MAX;
		while (i.good() && (l = *i).first.begin_ < k.end_) {
			if (l.second.count < cap)
				return 0;
			s = std::min(s, l.second.min_score);
			++i;
		}
		return s;
	}

	int max_score(interval k) const
	{
		Iterator i = begin(k.begin_);
		std::pair<interval, IntervalNode> l;
		int s = INT_MAX;
		while (i.good() && (l = *i).first.begin_ < k.end_) {
			s = std::min(s, l.second.max_score);
			++i;
		}
		assert(s != INT_MAX);
		return s;
	}

	Iterator begin(int p) const
	{
		const_iterator i = lower_bound(p), j;
		if (i == end() || i->first != p) {
			//assert(i != std::map<int, IntervalNode>::begin());
			j = i;
			i--;
		}
		else {
			j = i;
			++j;
		}
		return Iterator(i, j, *this);
	}

	const int cap;

};

#endif
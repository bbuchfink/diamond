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
#include <limits>
#include "interval.h"

struct IntervalPartition : protected std::map<int, int>
{

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
		std::pair<interval, int> operator*() const
		{
			return std::make_pair(interval(i_->first, j_ == parent_.end() ? std::numeric_limits<int>::max() : j_->first), i_->second);
		}
	private:
		const_iterator i_, j_;
		const IntervalPartition &parent_;
	};

	IntervalPartition()
	{
		(*this)[0] = 0;
	}

	void insert(interval k)
	{
		iterator i = lower_bound(k.begin_);
		if (i->first != k.begin_) {
			assert(i != begin());
			i--;
			i = std::map<int,int>::insert(std::make_pair(k.begin_, i->second)).first;
		}
		int last;
		while (i != end() && i->first < k.end_) {
			last = i->second++;
			++i;
		}
		if (i->first != k.end_)
			(*this)[k.end_] = last;
	}

	int covered(interval k, int n) const
	{
		Iterator i = begin(k.begin_);
		std::pair<interval, int> l;
		int c = 0;
		while (i.good() && (l = *i).first.begin_ < k.end_) {
			if (l.second >= n)
				c += k.overlap(l.first);
			++i;
		}
		return c;
	}

	Iterator begin(int p) const
	{
		const_iterator i = lower_bound(p), j;
		if (i->first != p) {
			assert(i != begin());
			j = i;
			i--;
		}
		else {
			j = i;
			++j;
		}
		return Iterator(i, j, *this);
	}

};

#endif
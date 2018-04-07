/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

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

#ifndef JOIN_RESULT_H_
#define JOIN_RESULT_H_

#include <algorithm>
#include <vector>
#include "../data_structures/double_array.h"

using std::pair;
using std::vector;

template<typename _t>
struct JoinResult : public vector<pair<DoubleArray<typename _t::Value>*, DoubleArray<typename _t::Value>*> >
{

	typedef typename _t::Value Value;

	struct Iterator
	{
		typename DoubleArray<Value>::Iterator r, s;
		Iterator(typename JoinResult<_t>::iterator begin, typename JoinResult<_t>::iterator end) :
			r(begin->first->begin()),
			s(begin->second->begin()),
			it_(begin),
			end_(end)
		{}
		Iterator& operator++()
		{
			++r;
			++s;
			if (!(r < it_->first->end())) {
				++it_;
				if (it_ < end_) {
					r = it_->first->begin();
					s = it_->second->begin();
				}
			}
			return *this;
		}
		bool good() const
		{
			return it_ < end_;
		}
	private:
		typename JoinResult<_t>::iterator it_, end_;
	};

	Iterator begin()
	{
		return Iterator(vector<pair<DoubleArray<typename _t::Value>*, DoubleArray<typename _t::Value>*> >::begin(), this->end());
	}

};

#endif
/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

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
#include <algorithm>
#include "../range.h"
#include "../data_structures/double_array.h"

template<typename _t>
struct JoinArrayIterator {

	JoinArrayIterator(_t *ptr, _t *end):
		ptr_(ptr),
		end_(end)
	{}

	Range<_t*> operator*() {
		return { ptr_ + 1, ptr_ + 1 + *ptr_ };
	}

	Range<_t*>* operator->() {
		range_ = this->operator*();
		return &range_;
	}

	JoinArrayIterator& operator++() {
		ptr_ += *ptr_ + 1;
		while (ptr_ < end_ && *ptr_ == 0)
			ptr_ += ptr_[1] + 1;
		range_ = { ptr_ + 1, ptr_ + 1 + *ptr_ };
		return *this;
	}
	
	operator bool() const {
		return ptr_ < end_;
	}

	void erase() {
		ptr_[1] = *ptr_;
		*ptr_ = 0;
	}

private:

	Range<_t*> range_;
	_t *ptr_, *end_;

};

template<typename _t>
struct JoinIterator {

	typename DoubleArray<_t>::Iterator r, s;

	JoinIterator(typename DoubleArray<_t>::Iterator r, typename DoubleArray<_t>::Iterator s):
		r(r),
		s(s)
	{}

	JoinIterator& operator++() {
		++r;
		++s;
		return *this;
	}

	operator bool() const {
		return r;
	}

	void erase() {
		r.erase();
		s.erase();
	}

};

/*template<typename _t>
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
*/

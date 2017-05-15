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

#ifndef GROWING_BUFFER_H_
#define GROWING_BUFFER_H_

#include <vector>

using std::vector;
using std::pair;

template<typename _t>
struct Growing_buffer
{

	inline void init(size_t size, size_t padding, size_t padding_front, _t init)
	{
		const size_t total = size + padding;
		data_.clear();
		data_.resize(total);
		col_size_ = total;
		for(size_t i=0;i<total;++i)
			data_[i] = init;
		init_ = init;
		center_.clear();
		center_.push_back(-1);
	}

	inline pair<_t*,_t*> get(int center)
	{
		data_.resize(data_.size() + col_size_);
		_t* ptr = last();
		for(size_t i=0;i<col_size_;++i)
			ptr[i] = init_;
		center_.push_back(center);
		return pair<_t*,_t*> (ptr-col_size_, ptr);
	}

	inline _t* last()
	{ return &*(data_.end() - col_size_); }

	const _t* column(int col) const
	{ return &data_[col_size_*col]; }

	int center(int col) const
	{ return center_[col]; }

private:
	vector<_t> data_;
	vector<int> center_;
	size_t col_size_;
	_t init_;

};

#endif

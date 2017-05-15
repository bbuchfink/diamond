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

#ifndef DOUBLE_BUFFER_H_
#define DOUBLE_BUFFER_H_

#include <vector>

using std::vector;
using std::pair;

template<typename _t>
struct Double_buffer
{

	inline void init(size_t size, size_t padding, size_t padding_front, _t init)
	{
		const size_t total = size + padding + padding_front;
		data_.clear();
		data_.resize(total*2);
		ptr1 = &data_[padding_front];
		ptr2 = &data_[total+padding_front];
		for(size_t i=0;i<total*2;++i)
			data_[i] = init;
	}

	inline pair<_t*,_t*> get(int)
	{ std::swap(ptr1, ptr2); return pair<_t*,_t*> (ptr2, ptr1); }

	inline _t* last()
	{ return ptr1; }

private:
	_t *ptr1, *ptr2;
	vector<_t> data_;

};

#endif /* DOUBLE_BUFFER_H_ */

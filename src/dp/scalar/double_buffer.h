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

template<typename T>
struct DoubleBuffer
{

	inline void init(size_t size, size_t padding, size_t padding_front, T init)
	{
		const size_t total = size + padding + padding_front;
		data_.clear();
		data_.resize(total * 2);
		ptr1 = &data_[padding_front];
		ptr2 = &data_[total + padding_front];
		for (size_t i = 0; i < total * 2; ++i)
			data_[i] = init;
	}

	inline std::pair<T*, T*> get(int)
	{
		std::swap(ptr1, ptr2); return std::pair<T*, T*>(ptr2, ptr1);
	}

	inline T* last()
	{
		return ptr1;
	}

private:
	T* ptr1, * ptr2;
	std::vector<T> data_;

};

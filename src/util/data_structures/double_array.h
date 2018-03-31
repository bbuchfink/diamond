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

#ifndef DOUBLE_ARRAY_H_
#define DOUBLE_ARRAY_H_

#include "../memory/memory_pool.h"

template<typename _t>
struct DoubleArray
{

	DoubleArray(unsigned n) :
		limits_(MemoryPool::global().alloc<unsigned>(n)),
		data_(NULL)
	{
	}

	~DoubleArray()
	{
		MemoryPool::global().free(limits_);
		MemoryPool::global().free(data_);
	}

	void init(unsigned data_size)
	{
		data_ = MemoryPool::global().alloc<_t>(data_size);
	}

	unsigned* limits()
	{
		return limits_;
	}

	_t* data()
	{
		return data_;
	}

private:

	unsigned *limits_;
	_t *data_;

};

#endif
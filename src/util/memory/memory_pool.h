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

#ifndef MEMORY_POOL_H_
#define MEMORY_POOL_H_

#include <list>
#include <map>
#include "../tinythread.h"

using std::list;
using std::map;

struct MemoryPool
{
	
	static void init(size_t size);
	static void* alloc(size_t n);
	static void free(void *p);

	template<typename _t>
	static _t* alloc(size_t n)
	{
		return (_t*)alloc(sizeof(_t)*n);
	}

private:

	struct Block
	{
		Block(size_t begin, size_t size) :
			begin(begin),
			size(size)
		{}
		size_t end() const
		{
			return begin + size;
		}
		size_t begin, size;
	};

	static list<Block>::iterator find_block(size_t min_size);

	static char *mem_;
	static list<Block> free_;
	static map<size_t, size_t> size_;
	static tthread::mutex mtx_;

};

#endif
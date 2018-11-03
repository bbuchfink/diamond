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

#ifndef MEMORY_POOL_H_
#define MEMORY_POOL_H_

#include <vector>
#include <algorithm>
#include <map>
#include "../tinythread.h"

using std::vector;
using std::pair;
using std::map;

struct Arena;

struct MemoryPool
{
	
	void* alloc(size_t n);
	void free(void *p);
	void clear();
	MemoryPool(bool thread_safe = true, size_t expected_limit = 0);
	void init(size_t expected_limit);

	size_t max_alloc_size()
	{
		return max_alloc_size_;
	}

	template<typename _t>
	_t* alloc(size_t n)
	{
		return (_t*)alloc(sizeof(_t)*n);
	}

	static MemoryPool& global()
	{
		return global_;
	}

private:

	enum { ARENA_SIZE_MULTIPLIER = 20 };

	typedef map<void*, pair<size_t, size_t> > SizeMap;

	tthread::mutex mtx_;
	vector<Arena*> arena_;
	SizeMap size_;
	size_t max_alloc_size_, current_alloc_size_, arena_size_;
	bool thread_safe_;

	static MemoryPool global_;

};

#endif
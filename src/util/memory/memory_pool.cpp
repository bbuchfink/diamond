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

#include <malloc.h>
#include <stdexcept>
#include "memory_pool.h"

char* MemoryPool::mem_ = NULL;
list<MemoryPool::Block> MemoryPool::free_;
map<size_t, size_t> MemoryPool::size_;
tthread::mutex MemoryPool::mtx_;

void MemoryPool::init(size_t size)
{
	mem_ = (char*)malloc(size);
	free_.push_back(Block(0, size));
}

list<MemoryPool::Block>::iterator MemoryPool::find_block(size_t min_size)
{
	list<Block>::iterator i;
	for (i = free_.begin(); i != free_.end(); ++i)
		if (i->size >= min_size)
			return i;
	return i;
}

void* MemoryPool::alloc(size_t n)
{
	mtx_.lock();
	list<Block>::iterator i = find_block(n);
	if (i == free_.end())
		throw std::runtime_error("MemoryPool: Failed to allocate memory.");
	void *m = mem_ + i->begin;
	size_[i->begin] = n;
	if (n == i->size)
		free_.erase(i);
	else {
		i->size -= n;
		i->begin += n;
	}
	mtx_.unlock();
	return m;
}

void MemoryPool::free(void *p)
{
	mtx_.lock();
	const size_t offset = (char*)p - mem_;
	map<size_t, size_t>::iterator s = size_.find(offset);
	if (s == size_.end())
		throw std::runtime_error("MemoryPool: Invalid free.");

	const size_t end = offset + s->second;
	list<Block>::iterator i;
	for (i = free_.begin(); i != free_.end(); ++i)
		if (i->begin >= end)
			break;

	list<Block>::iterator j = i;
	--j;
	if (i != free_.begin() && j->end() == offset) {
		j->size += s->second;
		if (i != free_.end() && i->begin == j->end()) {
			j->size += i->size;
			free_.erase(i);
		}
	}
	else if (i != free_.end() && end == i->begin) {
		i->begin -= s->second;
		i->size += s->second;
	} else
		free_.insert(i, Block(offset, s->second));

	size_.erase(s);
	mtx_.unlock();
}
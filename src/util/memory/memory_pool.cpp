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

#include <list>
#include <malloc.h>
#include <stdexcept>
#include "memory_pool.h"

using std::list;

vector<Arena*> MemoryPool::arena_;
tthread::mutex MemoryPool::mtx_;
MemoryPool::SizeMap MemoryPool::size_;
size_t MemoryPool::max_alloc_size_ = 0;
size_t MemoryPool::current_alloc_size_ = 0;
size_t MemoryPool::arena_size_ = 0;

struct Arena
{

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

	Arena(size_t size) :
		mem_((char*)malloc(size))
	{
		if (mem_ == NULL)
			throw std::runtime_error("Memory allocation failed due to insufficient memory.");
		free_.push_back(Block(0, size));
	}

	~Arena()
	{
		::free(mem_);
	}

	list<Block>::iterator find_block(size_t min_size)
	{
		list<Block>::iterator i;
		for (i = free_.begin(); i != free_.end(); ++i)
			if (i->size >= min_size)
				return i;
		return i;
	}

	void* alloc(size_t n)
	{
		list<Block>::iterator i = find_block(n);
		if (i == free_.end())
			return NULL;
		void *m = mem_ + i->begin;
		if (n == i->size)
			free_.erase(i);
		else {
			i->size -= n;
			i->begin += n;
		}
		return m;
	}

	void free(void *p, size_t size)
	{
		const size_t offset = (char*)p - mem_,
			end = offset + size;
		list<Block>::iterator i;
		for (i = free_.begin(); i != free_.end(); ++i)
			if (i->begin >= end)
				break;

		list<Block>::iterator j = i;
		--j;
		if (i != free_.begin() && j->end() == offset) {
			j->size += size;
			if (i != free_.end() && i->begin == j->end()) {
				j->size += i->size;
				free_.erase(i);
			}
		}
		else if (i != free_.end() && end == i->begin) {
			i->begin -= size;
			i->size += size;
		}
		else
			free_.insert(i, Block(offset, size));
	}

	char *mem_;
	list<Block> free_;

};

void MemoryPool::init(size_t expected_limit)
{
	arena_size_ = expected_limit / 20;
}

void* MemoryPool::alloc(size_t n)
{
	mtx_.lock();
	void *p;
	for (vector<Arena*>::iterator i = arena_.begin(); i < arena_.end(); ++i)
		if (p = (*i)->alloc(n)) {
			size_[p] = std::make_pair(i - arena_.begin(), n);
			mtx_.unlock();
			return p;
		}
	//const size_t alloc_size = n*ARENA_SIZE_MULTIPLIER;
	const size_t alloc_size = arena_size_;
	current_alloc_size_ += alloc_size;
	max_alloc_size_ = std::max(max_alloc_size_, current_alloc_size_);
	arena_.push_back(new Arena(alloc_size));
	p = arena_.back()->alloc(n);
	size_[p] = std::make_pair(arena_.size() - 1, n);
	mtx_.unlock();
	return p;
}

void MemoryPool::free(void *p)
{
	mtx_.lock();
	typename SizeMap::iterator s = size_.find(p);
	if (s == size_.end())
		throw std::runtime_error("MemoryPool: Invalid free.");
	arena_[s->second.first]->free(p, s->second.second);
	size_.erase(s);
	mtx_.unlock();
}

void MemoryPool::clear()
{
	for (vector<Arena*>::iterator i = arena_.begin(); i < arena_.end(); ++i)
		delete *i;
	arena_.clear();
	current_alloc_size_ = 0;
}
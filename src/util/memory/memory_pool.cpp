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
#include <stdlib.h>
#include <stdexcept>
#include "memory_pool.h"

using std::list;

MemoryPool MemoryPool::global_;

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

	void* alloc(size_t n)
	{
		if (free_.empty())
			return NULL;
		Block &b = free_.back();
		if (b.size < n)
			return NULL;
		void *m = mem_ + b.begin;
		if (n == b.size)
			free_.pop_back();
		else {
			b.size -= n;
			b.begin += n;
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

MemoryPool::MemoryPool(bool thread_safe, size_t expected_limit):
	max_alloc_size_(0),
	current_alloc_size_(0),
	arena_size_(expected_limit ? expected_limit / 100 : 256 * (1<<20)),
	thread_safe_(thread_safe)
{
}

void MemoryPool::init(size_t expected_limit)
{

}


void* MemoryPool::alloc(size_t n)
{
	if (thread_safe_) mtx_.lock();
	void *p;
	for (vector<Arena*>::iterator i = arena_.begin(); i < arena_.end(); ++i)
		if ((p = (*i)->alloc(n))) {
			size_[p] = std::make_pair(i - arena_.begin(), n);
			if (thread_safe_) mtx_.unlock();
			return p;
		}
	const size_t alloc_size = std::max(arena_size_, n*ARENA_SIZE_MULTIPLIER);
	current_alloc_size_ += alloc_size;
	max_alloc_size_ = std::max(max_alloc_size_, current_alloc_size_);
	arena_.push_back(new Arena(alloc_size));
	p = arena_.back()->alloc(n);
	size_[p] = std::make_pair(arena_.size() - 1, n);
	if (thread_safe_) mtx_.unlock();
	return p;
}

void MemoryPool::free(void *p)
{
	if(thread_safe_) mtx_.lock();
	SizeMap::iterator s = size_.find(p);
	if (s == size_.end())
		throw std::runtime_error("MemoryPool: Invalid free.");
	arena_[s->second.first]->free(p, s->second.second);
	size_.erase(s);
	if (thread_safe_) mtx_.unlock();
}

void MemoryPool::clear()
{
	for (vector<Arena*>::iterator i = arena_.begin(); i < arena_.end(); ++i)
		delete *i;
	arena_.clear();
	current_alloc_size_ = 0;
}
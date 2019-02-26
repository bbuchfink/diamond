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

#include <limits>
#include <mutex>

struct Queue
{
	enum { end = size_t(-1) };
	Queue(size_t begin, size_t end) :
		next_(begin),
		block_(false),
		end_(end)
	{}
	template<typename _f>
	size_t get(_f &f)
	{
		std::unique_lock<std::mutex> lock(mtx_);
		while (block_)
			cond_.wait(lock);
		const size_t q = next_++;
		if (q >= end_) {
			return Queue::end;
		}
		block_ = f(q);
		return q;
	}
	size_t next() const
	{
		return next_;
	}
	size_t get_end() const
	{
		return end_;
	}
	void release() {
		block_ = false;
		cond_.notify_all();
	}
private:
	std::mutex mtx_;
	std::condition_variable cond_;
	volatile size_t next_;
	volatile bool block_;
	const size_t end_;
};
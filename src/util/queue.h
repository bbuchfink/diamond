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

#include <limits>
#include "tinythread.h"

struct Queue
{
	enum { end = std::numeric_limits<size_t>::max() };
	Queue(size_t begin, size_t end) :
		next_(begin),
		end_(end)
	{}
	template<typename _f>
	size_t get(_f &f)
	{
		mtx_.lock();
		const size_t q = next_++;
		if (q >= end_) {
			mtx_.unlock();
			return Queue::end;
		}
		f(q);
		mtx_.unlock();
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
private:
	tthread::mutex mtx_;
	size_t next_;
	const size_t end_;
};
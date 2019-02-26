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

#ifndef THREAD_H_
#define THREAD_H_

#include <vector>
#include <exception>
#include <stdexcept>
#include <mutex>

using std::vector;

template<typename _t>
struct Atomic
{
	Atomic(const _t &v):
		v_ (v)
	{ }
	Atomic& operator=(const _t &v)
	{
		v_ = v;
		return *this;
	}
	volatile _t operator++(int)
	{
		mtx_.lock();
		_t r = v_++;
		mtx_.unlock();
		return r;
	}
	_t operator--(int)
	{
		mtx_.lock();
		_t r = v_--;
		mtx_.unlock();
		return r;
	}
	_t post_add(const _t &v) {
		mtx_.lock();
		_t r = v_;
		v_ += v;
		mtx_.unlock();
		return r;
	}
private:
	volatile _t v_;
	std::mutex mtx_;
};

#endif /* THREAD_H_ */

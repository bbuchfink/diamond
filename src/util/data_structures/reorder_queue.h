/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once
#include <mutex>

template<typename T, typename F>
struct ReorderQueue
{
	ReorderQueue(size_t begin, F& f) :
		f_(f),
		begin_(begin),
		next_(begin),
		size_(0),
		max_size_(0)
	{}

	size_t size() const
	{
		return size_;
	}
	size_t max_size() const
	{
		return max_size_;
	}
	size_t next() const
	{
		return next_;
	}
	size_t begin() const {
		return begin_;
	}

	void push(size_t n, T value)
	{
		mtx_.lock();
		//cout << "n=" << n << " next=" << next_ << endl;
		if (n != next_) {
			backlog_[n] = value;
			size_ += value ? value->alloc_size() : 0;
			max_size_ = std::max(max_size_, size_);
			mtx_.unlock();
		}
		else
			flush(value);
	}

private:

	void flush(T value)
	{
		size_t n = next_ + 1;
		std::vector<T> out;
		out.push_back(value);
		typename std::map<size_t, T>::iterator i;
		do {
			while ((i = backlog_.begin()) != backlog_.end() && i->first == n) {
				out.push_back(i->second);
				backlog_.erase(i);
				++n;
			}
			mtx_.unlock();
			size_t size = 0;
			for (typename std::vector<T>::iterator j = out.begin(); j < out.end(); ++j) {
				if (*j) {
					f_(*j);
					if (*j != value)
						size += (*j)->alloc_size();
					delete *j;
				}
			}
			out.clear();
			mtx_.lock();
			size_ -= size;
		} while ((i = backlog_.begin()) != backlog_.end() && i->first == n);
		next_ = n;
		mtx_.unlock();
	}

	std::mutex mtx_;
	F& f_;
	std::map<size_t, T> backlog_;
	size_t begin_, next_, size_, max_size_;
};
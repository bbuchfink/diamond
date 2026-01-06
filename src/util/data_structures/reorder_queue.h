/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

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
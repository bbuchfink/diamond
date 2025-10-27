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

#include <list>
#include <vector>
#include <iterator>
#include <stddef.h>
#include "../parallel/mutex.h"

template<typename T, size_t E, typename Sync = ::Sync>
struct Deque {

	static const size_t EXPONENT = E;
	static const size_t N = (size_t)1 << E;
	static const ptrdiff_t SHIFT = E;
	static const ptrdiff_t MASK = ((ptrdiff_t)1 << E) - 1;

	typedef std::vector<T> Bucket;

	Deque()
	{
		new_bucket();
	}

	void reserve(size_t n) {}

	void push_back(const T& x) {
		if (buckets.back().size() >= N) {
			new_bucket();
		}
		buckets.back().push_back(x);
	}

	void push_back(const T* ptr, size_t n) {
		if (buckets.back().size() + n > N) {
			new_bucket();
		}
		buckets.back().insert(buckets.back().end(), ptr, ptr + n);
	}

	template<typename It>
	void push_back(It begin, It end) {
		mtx_.lock();
		while (begin < end) {
			if (buckets.back().size() == N)
				new_bucket();
			const ptrdiff_t n = std::min(end - begin, ptrdiff_t(N - buckets.back().size()));
			buckets.back().insert(buckets.back().end(), begin, begin + n);
			begin += n;
		}		
		mtx_.unlock();
	}

	size_t size() const {
		size_t n = 0;
		for (const Bucket& b : buckets)
			n += b.size();
		return n;
	}

	void move(std::vector<T>& dst) {
		if (buckets.size() == 1 && dst.empty())
			dst = std::move(buckets.front());
		else {
			for (const Bucket& b : buckets)
				dst.insert(dst.end(), b.begin(), b.end());
		}
		buckets.clear();
	}

	struct Iterator {

		using iterator_category = std::random_access_iterator_tag;
		using difference_type = ptrdiff_t;
		using value_type = T;
		using pointer = T*;
		using reference = T&;

		Iterator() {}

		Iterator(size_t i, T** data) :
			i_(i),
			data_(data)
		{
		}

		T& operator*() {
			return data_[i_ >> SHIFT][i_ & MASK];
		}

		T* operator->() {
			return &data_[i_ >> SHIFT][i_ & MASK];
		}

		T& operator*() const {
			return data_[i_ >> SHIFT][i_ & MASK];
		}

		T& operator[](ptrdiff_t i) {
			ptrdiff_t j = i_ + i;
			return data_[j >> SHIFT][j & MASK];
		}

		ptrdiff_t operator-(Iterator& it) {
			return i_ - it.i_;
		}

		ptrdiff_t operator-(const Iterator& it) const {
			return i_ - it.i_;
		}

		Iterator operator+(ptrdiff_t i) const {
			return Iterator(i_ + i, data_);
		}

		Iterator operator-(ptrdiff_t i) const {
			return Iterator(i_ - i, data_);
		}

		Iterator& operator++() {
			++i_;
			return *this;
		}

		Iterator operator++(int) {
			Iterator r(i_, data_);
			++i_;
			return r;
		}

		Iterator& operator--() {
			--i_;
			return *this;
		}

		bool operator==(const Iterator& it) const {
			return i_ == it.i_;
		}

		bool operator!=(const Iterator& it) const {
			return i_ != it.i_;
		}

		bool operator>=(const Iterator& it) const {
			return i_ >= it.i_;
		}

		bool operator<=(const Iterator& it) const {
			return i_ <= it.i_;
		}

		bool operator<(const Iterator& it) const {
			return i_ < it.i_;
		}

		bool operator>(const Iterator& it) const {
			return i_ > it.i_;
		}

#if SIZEOF_INT != SIZEOF_PTRDIFF_T
		Iterator operator-(int i) const {
			return Iterator(i_ - i, data_);
		}
#endif

		Iterator& operator+=(ptrdiff_t i) {
			i_ += i;
			return *this;
		}

		Iterator& operator-=(ptrdiff_t i) {
			i_ -= i;
			return *this;
		}

	private:

		ptrdiff_t i_;
		T** data_;

	};

	Iterator begin() {
		init();
		return Iterator(0, data_.data());
	}

	Iterator end() {
		init();
		return Iterator(total_, data_.data());
	}

private:

	void init() {
		data_.clear();
		total_ = 0;
		for (Bucket& b : buckets) {
			data_.push_back(b.data());
			total_ += b.size();
		}
	}

	void new_bucket() {
		buckets.emplace_back();
		buckets.back().reserve(N);
	}

	std::list<Bucket> buckets;
	std::vector<T*> data_;
	size_t total_;
	Mutex<Sync> mtx_;

};

template<typename T, size_t E>
struct AsyncWriter {

	AsyncWriter(Deque<T, E, Async>& dst):
		dst_(&dst)
	{}

	void write(const T& v) {
		buf_.push_back(v);
		if (buf_.size() >= BUF_SIZE) {
			dst_->push_back(buf_.begin(), buf_.end());
			buf_.clear();
		}
	}

	~AsyncWriter() {
		dst_->push_back(buf_.begin(), buf_.end());
	}

private:

	static const size_t BUF_SIZE = 4096;

	Deque<T, E, Async>* dst_;
	std::vector<T> buf_;

};
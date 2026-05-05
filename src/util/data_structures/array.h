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
#include <stdint.h>
#include <algorithm>

template<typename T>
struct Array {

	using Size = int64_t;

	Array():
		ptr_(nullptr),
		size_(0)
	{}

	Array(Size alloc_size) :
		ptr_(new T[alloc_size]),
		size_(0)
	{}

	Array& operator=(Array&& a) noexcept {
		ptr_ = a.ptr_;
		size_ = a.size_;
		a.ptr_ = nullptr;
		return *this;
	}

	~Array() {
		delete[] ptr_;
	}

	T* data() {
		return ptr_;
	}

	T* begin() {
		return ptr_;
	}

	T* end() {
		return ptr_ + size_;
	}

	Size size() const {
		return size_;
	}

	void assign(const T& x) {
		*ptr_ = x;
		size_ = 1;
	}

	template<typename It>
	void assign(It begin, It end) {
		std::copy(begin, end, ptr_);
		size_ = end - begin;
	}

	template<typename It>
	void assign_reversed(It begin, It end) {
		T* p = ptr_;
		for (It i = end - 1; i >= begin; --i)
			*p++ = *i;
		size_ = end - begin;
	}

	template<typename It>
	void push_back(It begin, It end) {
		std::copy(begin, end, this->end());
		size_ += end - begin;
	}

	template<typename It>
	void push_back_reversed(It begin, It end) {
		T* p = this->end();
		for (It i = end - 1; i >= begin; --i)
			*p++ = *i;
		size_ += end - begin;
	}

	void push_back(Size n, const T& value) {
		T* p = end();
		std::fill(p, p + n, value);
		size_ += n;
	}

	T operator[](int i) const {
		return ptr_[i];
	}

private:
	T* ptr_;
	Size size_;

};
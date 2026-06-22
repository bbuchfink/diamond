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
#include <stddef.h>
#include <stdint.h>
#include <algorithm>
#include "../range.h"

template<typename T>
struct DoubleArray {

	static constexpr size_t header_size = sizeof(uint32_t);

	DoubleArray()
	{}

	DoubleArray(void* data, size_t size = 0):
		data_((char*)data),
		size_(size)
	{}

	struct Iterator {

		Iterator(void *ptr) :
			ptr_((char*)ptr),
			end_(nullptr)
		{}

		Iterator(char *ptr, char *end) :
			ptr_(ptr),
			end_(end)
		{
			skip_del();
		}

		uint32_t& count() {
			return *(uint32_t*)ptr_;
		}

		Range<T*> operator*() {
			return { (T*)(ptr_ + header_size), (T*)(ptr_ + header_size) + count() };
		}

		Range<T*>* operator->() {
			range_ = this->operator*();
			return &range_;
		}

		Iterator& operator++() {
			next();
			skip_del();
			return *this;
		}

		void next() {
			ptr_ += count() * sizeof(T) + header_size;
		}

		operator bool() const {
			return ptr_ < end_;
		}

		void erase() {
			uint32_t n = count();
			*(uint32_t*)(ptr_ + header_size) = n;
			count() = 0;
			ptr_ += n * sizeof(T) + header_size;
		}

		ptrdiff_t operator-(const Iterator &x) const {
			return ptr_ - x.ptr_;
		}

		bool operator==(const Iterator &x) const {
			return ptr_ == x.ptr_;
		}

	private:

		void skip_del() {
			while (ptr_ < end_ && count() == 0)
				ptr_ += (*(uint32_t*)(ptr_ + header_size)) * sizeof(T) + header_size;
		}

		Range<T*> range_;
		char *ptr_, *end_;

		friend struct DoubleArray;

	};

	Iterator begin() {
		return Iterator(data_, data_ + size_);
	}

	void set_end(const Iterator &it) {
		size_ = it.ptr_ - data_;
	}

	void append(const DoubleArray &d) {
		if (d.data_ == data_ + size_) {
			size_ += d.size_;
			return;
		}
		std::copy(d.data_, d.data_ + d.size_, data_ + size_);
		size_ += d.size_;
	}

	size_t offset(const Iterator &it) const {
		return size_t(it.ptr_ - data_);
	}

	T& operator[](size_t i) {
		return *(T*)(data_ + i);
	}

private:

	char *data_;
	size_t size_;

};

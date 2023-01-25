/****
DIAMOND protein aligner
Copyright (C) 2013-2021 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen

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

#pragma once
#include <stdint.h>
#include <stddef.h>
#include <algorithm>
#include "../range.h"

template<typename T>
struct DoubleArray {

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
			return { (T*)(ptr_ + 4), (T*)(ptr_ + 4) + count() };
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
			ptr_ += count() * sizeof(T) + 4;
		}

		operator bool() const {
			return ptr_ < end_;
		}

		void erase() {
			uint32_t n = count();
			*(uint32_t*)(ptr_ + 4) = n;
			count() = 0;
			ptr_ += n * sizeof(T) + 4;
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
				ptr_ += (*(uint32_t*)(ptr_ + 4)) * sizeof(T) + 4;
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

	uint32_t offset(const Iterator &it) const {
		return uint32_t(it.ptr_ - data_);
	}

	T& operator[](uint32_t i) {
		return *(T*)(data_ + i);
	}

private:

	char *data_;
	size_t size_;

};

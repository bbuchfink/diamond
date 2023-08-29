/****
DIAMOND protein aligner
Copyright (C) 2020 Max Planck Society for the Advancement of Science e.V.

Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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
#include "../memory/alignment.h"

template<typename T>
struct MemBuffer {

	enum { ALIGN = 32 };

	typedef T value_type;

	MemBuffer():
		data_(nullptr),
		size_(0),
		alloc_size_(0)
	{}

	MemBuffer(size_t n):
		data_((T*)Util::Memory::aligned_malloc(n * sizeof(T), ALIGN)),
		size_(n),
		alloc_size_(n)
	{
	}

	~MemBuffer() {
		Util::Memory::aligned_free(data_);
	}

	void resize(size_t n) {
		if (alloc_size_ < n) {
			Util::Memory::aligned_free(data_);
			data_ = (T*)Util::Memory::aligned_malloc(n * sizeof(T), ALIGN);
			alloc_size_ = n;
		}
		size_ = n;
	}

	size_t size() const {
		return size_;
	}

	T* begin() {
		return data_;
	}

	T* end() {
		return data_ + size_;
	}

	const T* begin() const {
		return data_;
	}

	const T* end() const {
		return data_ + size_;
	}

	T& operator[](size_t i) {
		return data_[i];
	}

	const T& operator[](size_t i) const {
		return data_[i];
	}

private:

	T *data_;
	size_t size_, alloc_size_;

};
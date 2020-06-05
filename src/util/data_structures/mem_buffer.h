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
#include <stdlib.h>
#include <exception>
#include "../memory/alignment.h"

template<typename _t>
struct MemBuffer {

	enum { ALIGN = 32 };

	MemBuffer():
		data_(nullptr),
		size_(0),
		alloc_size_(0)
	{}

	MemBuffer(size_t n):
		data_((_t*)aligned_malloc(n * sizeof(_t), ALIGN)),
		size_(n),
		alloc_size_(n)
	{
		if (data_ == nullptr)
			throw std::bad_alloc();
	}

	~MemBuffer() {
		aligned_free(data_);
	}

	void resize(size_t n) {
		if (alloc_size_ < n) {
			aligned_free(data_);
			data_ = (_t*)aligned_malloc(n * sizeof(_t), ALIGN);
			if (data_ == nullptr)
				throw std::bad_alloc();
			alloc_size_ = n;
		}
		size_ = n;
	}

	size_t size() const {
		return size_;
	}

	_t* begin() {
		return data_;
	}

	_t* end() {
		return data_ + size_;
	}

	_t& operator[](size_t i) {
		return data_[i];
	}

	const _t& operator[](size_t i) const {
		return data_[i];
	}

private:

	_t *data_;
	size_t size_, alloc_size_;

};
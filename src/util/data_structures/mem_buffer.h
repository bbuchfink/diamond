/****
DIAMOND protein aligner
Copyright (C) 2013-2019 Benjamin Buchfink <buchfink@gmail.com>

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


#ifndef MEM_BUFFER_H_
#define MEM_BUFFER_H_

#include <stdlib.h>

template<typename _t>
struct MemBuffer {

	MemBuffer():
		data_(nullptr),
		size_(0),
		alloc_size_(0)
	{}

	MemBuffer(size_t n):
		data_((_t*)malloc(n * sizeof(_t))),
		size_(n),
		alloc_size_(n)
	{}

	~MemBuffer() {
		free(data_);
	}

	void resize(size_t n) {
		if (alloc_size_ < n) {
			free(data_);
			data_ = (_t*)malloc(n * sizeof(_t));
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

private:

	_t *data_;
	size_t size_, alloc_size_;

};

#endif
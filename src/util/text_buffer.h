/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#ifndef TEXT_BUFFER_H_
#define TEXT_BUFFER_H_

#include <stdlib.h>
#include <stdio.h>
#include <stdexcept>
#include <stdint.h>
#include <limits>
#include "util.h"

struct Text_buffer
{

	Text_buffer():
		data_ (0),
		ptr_ (data_),
		alloc_size_(0)
	{ }

	~Text_buffer()
	{
		free(data_);
	}

	void reserve(size_t n)
	{
		const size_t s = ptr_ - data_;
		alloc_size_ = s + n + block_size - ((s + n) & (block_size - 1));
		data_ = (char*)realloc(data_, alloc_size_);
		ptr_ = data_ + s;
		if (data_ == 0) throw std::runtime_error("Failed to allocate memory.");
	}

	void operator+=(size_t n)
	{ ptr_ += n; }

	operator char*() const
	{ return ptr_; }

	char* get_begin() const
	{ return data_; }

	void clear()
	{ ptr_ = data_; }

	template<typename _t>
	Text_buffer& write(const _t& data)
	{
		reserve(sizeof(_t));
		*reinterpret_cast<_t*>(ptr_) = data;
		ptr_ += sizeof(_t);
		return *this;
	}

	void write_c_str(const char* s)
	{
		const size_t l = strlen(s)+1;
		reserve(l);
		memcpy(ptr_, s, l);
		ptr_ += l;
	}

	void write_c_str(const char*s, size_t len)
	{
		reserve(len+1);
		memcpy(ptr_, s, len);
		ptr_ += len;
		*(ptr_++) = '\0';
	}

	void write_raw(const char *s, size_t len)
	{
		reserve(len);
		memcpy(ptr_, s, len);
		ptr_ += len;
	}

	void write_until(const char *s, const char *delimiters)
	{
		write_raw(s, find_first_of(s, delimiters));
	}

	Text_buffer& write_packed(unsigned x)
	{
		if(x <= (unsigned)std::numeric_limits<uint8_t>::max())
			write((uint8_t)x);
		else if(x <= (unsigned)std::numeric_limits<uint16_t>::max())
			write((uint16_t)x);
		else
			write(x);
		return *this;
	}

	size_t size() const
	{ return ptr_ - data_; }

	size_t alloc_size() const
	{
		return alloc_size_;
	}

	Text_buffer& operator<<(const string &s)
	{
		const size_t l = s.length();
		reserve(l);
		memcpy(ptr_, s.c_str(), l);
		ptr_ += l;
		return *this;
	}

	Text_buffer& operator<<(const char* s)
	{
		const size_t l = strlen(s);
		reserve(l);
		memcpy(ptr_, s, l);
		ptr_ += l;
		return *this;
	}

	Text_buffer& operator<<(char c)
	{
		reserve(1);
		*(ptr_++) = c;
		return *this;
	}
	
	Text_buffer& operator<<(unsigned int x)
	{
		//write(x);
		reserve(16);
		ptr_ += sprintf(ptr_, "%u", x);
		return *this;
	}

	Text_buffer& operator<<(int x)
	{
		//write(x);
		reserve(16);
		ptr_ += sprintf(ptr_, "%i", x);
		return *this;
	}

	Text_buffer& operator<<(unsigned long x)
	{
		reserve(32);
		ptr_ += sprintf(ptr_, "%lu", x);
		return *this;
	}
	
	Text_buffer& operator<<(unsigned long long x)
	{
		reserve(32);
		ptr_ += sprintf(ptr_, "%llu", x);
		return *this;
	}

	Text_buffer& operator<<(double x)
	{
		reserve(32);
		ptr_ += sprintf(ptr_, "%.1lf", x);
		return *this;
	}

	Text_buffer& print_d(double x)
	{
		reserve(32);
		ptr_ += sprintf(ptr_, "%lf", x);
		return *this;
	}

	Text_buffer& print_e(double x)
	{
		reserve(32);
		ptr_ += sprintf(ptr_, "%.1le", x);
		return *this;
	}

	Text_buffer& print(unsigned i, unsigned width)
	{
		reserve(16);
		char buf[16];
		const int n = sprintf(buf, "%u", i),
			padding = std::max((int)width - n, 0);
		memset(ptr_, ' ', padding);
		memcpy(ptr_ + padding, buf, n);
		ptr_ += padding + n;
		return *this;
	}

	/*Text_buffer& operator<<(uint8_t x)
	{
		write(x);
		return *this;
	}*/

	template<typename _t>
	Text_buffer& operator<<(const vector<_t> &v)
	{
		const size_t l = v.size() * sizeof(_t);
		reserve(l);
		memcpy(ptr_, v.data(), l);
		ptr_ += l;
		return *this;
	}

	char& operator[](size_t pos)
	{
		return data_[pos];
	}

protected:
	enum { block_size = 4096 };
	char *data_, *ptr_;
	size_t alloc_size_;

};

#endif /* TEXT_BUFFER_H_ */

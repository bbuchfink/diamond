/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen
						
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
#include <stdio.h>
#include <stdexcept>
#include <stdint.h>
#include <limits>
#include <vector>
#include "algo/varint.h"
#include "string/string.h"

inline size_t find_first_of(const char* s, const char* delimiters)
{
	const char* t = s;
	while (*t && strchr(delimiters, *t) == 0)
		++t;
	return t - s;
}

struct TextBuffer
{

	typedef char value_type;

	TextBuffer():
		data_ (0),
		ptr_ (data_),
		alloc_size_(0)
	{
	}

	~TextBuffer()
	{
		free(data_);
	}

	void reserve(size_t n)
	{
		const size_t s = ptr_ - data_;
		if (s + n < alloc_size_)
			return;
		alloc_size_ = s + n + block_size - ((s + n) & (block_size - 1));
		data_ = (char*)realloc(data_, alloc_size_);
		ptr_ = data_ + s;
		if (data_ == 0) throw std::runtime_error("Failed to allocate memory.");
	}

	void operator+=(size_t n)
	{ ptr_ += n; }

	operator char*() const
	{ return ptr_; }

	char* data() const
	{ return data_; }

	void clear()
	{ ptr_ = data_; }

	template<typename _t>
	TextBuffer& write(const _t& data)
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

	void push_back(char c) {
		*(ptr_++) = c;
	}

	void write_until(const char *s, const char *delimiters)
	{
		write_raw(s, find_first_of(s, delimiters));
	}

	TextBuffer& write_packed(unsigned x)
	{
		if(x <= (unsigned)std::numeric_limits<uint8_t>::max())
			write((uint8_t)x);
		else if(x <= (unsigned)std::numeric_limits<uint16_t>::max())
			write((uint16_t)x);
		else
			write(x);
		return *this;
	}

	TextBuffer& write_varint(unsigned x)
	{
		::write_varint(x, *this);
		return *this;
	}

	size_t size() const
	{ return ptr_ - data_; }

	size_t alloc_size() const
	{
		return alloc_size_;
	}

	TextBuffer& operator<<(const std::string &s)
	{
		const size_t l = s.length();
		reserve(l);
		memcpy(ptr_, s.c_str(), l);
		ptr_ += l;
		return *this;
	}

	TextBuffer& operator<<(const char* s)
	{
		const size_t l = strlen(s);
		reserve(l);
		memcpy(ptr_, s, l);
		ptr_ += l;
		return *this;
	}

	TextBuffer& operator<<(char c)
	{
		reserve(1);
		*(ptr_++) = c;
		return *this;
	}

	TextBuffer& operator<<(unsigned char c)
	{
		reserve(1);
		*(ptr_++) = c;
		return *this;
	}
	
	TextBuffer& operator<<(unsigned int x)
	{
		//write(x);
		reserve(16);
		ptr_ += sprintf(ptr_, "%u", x);
		return *this;
	}

	TextBuffer& operator<<(int x)
	{
		//write(x);
		reserve(16);
		ptr_ += sprintf(ptr_, "%i", x);
		return *this;
	}

	TextBuffer& operator<<(unsigned long x)
	{
		reserve(32);
		ptr_ += sprintf(ptr_, "%lu", x);
		return *this;
	}
	
	TextBuffer& operator<<(unsigned long long x)
	{
		reserve(32);
		ptr_ += sprintf(ptr_, "%llu", x);
		return *this;
	}

	TextBuffer& operator<<(long x)
	{
		reserve(32);
		ptr_ += sprintf(ptr_, "%li", x);
		return *this;
	}

	TextBuffer& operator<<(long long x)
	{
		reserve(32);
		ptr_ += sprintf(ptr_, "%lli", x);
		return *this;
	}

	TextBuffer& operator<<(double x)
	{
		reserve(32);
		ptr_ += Util::String::format_double(x, ptr_);
		return *this;
	}

	TextBuffer& print_d(double x)
	{
		reserve(32);
		ptr_ += sprintf(ptr_, "%lf", x);
		return *this;
	}

	TextBuffer& print_e(double x)
	{
		reserve(32);
		if (x == 0.0)
			ptr_ += sprintf(ptr_, "0.0");
		else
			ptr_ += sprintf(ptr_, "%.2e", x);
		return *this;
	}

	TextBuffer& print(unsigned i, unsigned width)
	{
		reserve(16);
		sprintf(ptr_, "%*u", width, i);
		ptr_ += width;
		return *this;
	}

	template<typename T>
	TextBuffer& print(const std::vector<T> &v, char separator)
	{
		if (v.empty()) return *this;
		typename std::vector<T>::const_iterator i = v.begin();
		*this << *(i++);
		for (; i < v.end(); ++i)
			*this << ';' << *i;
		return *this;
	}

	template<typename T>
	TextBuffer& operator<<(const std::vector<T> &v)
	{
		const size_t l = v.size() * sizeof(T);
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

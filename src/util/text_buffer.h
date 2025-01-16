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
#include <stdio.h>
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

	using value_type = char;

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
		alloc_size_ = s + n + BLOCK_SIZE - ((s + n) & (BLOCK_SIZE - 1));
		data_ = (char*)realloc(data_, alloc_size_);
		if (data_ == nullptr) throw std::bad_alloc();
		ptr_ = data_ + s;
	}

	void operator+=(size_t n)
	{ ptr_ += n; }

	operator char*() const
	{ return ptr_; }

	char* data() const
	{ return data_; }

	void clear()
	{ ptr_ = data_; }

	template<typename T>
	TextBuffer& write(const T& data)
	{
		reserve(sizeof(T));
#ifdef __sparc__
		std::copy((char*)&data, (char*)&data + sizeof(T), ptr_);
#else
		*reinterpret_cast<T*>(ptr_) = data;
#endif
		ptr_ += sizeof(T);
		return *this;
	}

	void write_c_str(const char* s)
	{
		const size_t l = strlen(s) + 1;
		reserve(l);
		memcpy(ptr_, s, l);
		ptr_ += l;
	}

	void write_c_str(const char*s, size_t len)
	{
		reserve(len + 1);
		memcpy(ptr_, s, len);
		ptr_ += len;
		*ptr_++ = '\0';
	}

	void write_raw(const char *s, size_t len)
	{
		reserve(len);
		memcpy(ptr_, s, len);
		ptr_ += len;
	}

	void push_back(char c) {
		*ptr_++ = c;
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
		char buf[5];
		char* end = write_varuint32(x, buf);
		write_raw(buf, end - buf);
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
		*ptr_++ = c;
		return *this;
	}

	TextBuffer& operator<<(unsigned char c)
	{
		reserve(1);
		*ptr_++ = c;
		return *this;
	}
	
	TextBuffer& operator<<(unsigned int x)
	{
		//write(x);
		reserve(16);
		ptr_ += snprintf(ptr_, 16, "%u", x);
		return *this;
	}

	TextBuffer& operator<<(int x)
	{
		//write(x);
		reserve(16);
		ptr_ += snprintf(ptr_, 16, "%i", x);
		return *this;
	}

	TextBuffer& operator<<(unsigned long x)
	{
		reserve(32);
		ptr_ += snprintf(ptr_, 32, "%lu", x);
		return *this;
	}
	
	TextBuffer& operator<<(unsigned long long x)
	{
		reserve(32);
		ptr_ += snprintf(ptr_, 32, "%llu", x);
		return *this;
	}

	TextBuffer& operator<<(long x)
	{
		reserve(32);
		ptr_ += snprintf(ptr_, 32, "%li", x);
		return *this;
	}

	TextBuffer& operator<<(long long x)
	{
		reserve(32);
		ptr_ += snprintf(ptr_, 32, "%lli", x);
		return *this;
	}

	TextBuffer& operator<<(double x)
	{
		reserve(32);
		ptr_ += Util::String::format_double(x, ptr_, 32);
		return *this;
	}

	TextBuffer& print_d(double x)
	{
		reserve(32);
		ptr_ += snprintf(ptr_, 32, "%lf", x);
		return *this;
	}

	TextBuffer& print_e(double x)
	{
		reserve(32);
		if (x == 0.0)
			ptr_ += snprintf(ptr_, 32, "0.0");
		else
			ptr_ += snprintf(ptr_, 32, "%.2e", x);
		return *this;
	}

	TextBuffer& print(unsigned i, unsigned width)
	{
		reserve(16);
		snprintf(ptr_, 16, "%*u", width, i);
		ptr_ += width;
		return *this;
	}

	template<typename T>
	TextBuffer& print(const std::vector<T> &v, char separator)
	{
		if (v.empty()) return *this;
		typename std::vector<T>::const_iterator i = v.begin();
		*this << *i++;
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
	enum { BLOCK_SIZE = 4096 };
	char *data_, *ptr_;
	size_t alloc_size_;

};

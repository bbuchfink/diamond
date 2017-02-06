/****
Copyright (c) 2014-2016, University of Tuebingen, Benjamin Buchfink
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
		ptr_ (data_)
	{ }

	~Text_buffer()
	{
		free(data_);
	}

	void reserve(size_t n)
	{
		const size_t s = ptr_ - data_, new_size = s + n + block_size - ((s + n) & (block_size - 1));
		data_ = (char*)realloc(data_, new_size);
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
		ptr_ += sprintf(ptr_, "%4u", i);
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
	enum { block_size = 65536 };
	char *data_, *ptr_;

};

#endif /* TEXT_BUFFER_H_ */

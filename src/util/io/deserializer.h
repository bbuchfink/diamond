/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <string.h>
#include <utility>
#include <iterator>
#include <assert.h>
#include "input_stream_buffer.h"
#include "../system/endianness.h"

struct Deserializer {

	Deserializer(InputStreamBuffer* buffer);
	void rewind();
	Deserializer& seek(int64_t pos);
	void seek_forward(size_t n);
	bool seek_forward(char delimiter);
	void close();

	Deserializer& operator>>(unsigned &x)
	{
		read(x);
		x = big_endian_byteswap(x);
		return *this;
	}

	Deserializer& operator>>(int32_t &x)
	{
		read(x);
		x = big_endian_byteswap(x);
		return *this;
	}

	Deserializer& operator>>(int64_t& x)
	{
		read(x);
		x = big_endian_byteswap(x);
		return *this;
	}

	Deserializer& operator>>(short& x) {
		read(x);
		x = big_endian_byteswap(x);
		return *this;
	}

	Deserializer& operator>>(unsigned short& x) {
		read(x);
		x = big_endian_byteswap(x);
		return *this;
	}

	Deserializer& operator>>(unsigned long &x)
	{
		read(x);
		x = big_endian_byteswap(x);
		return *this;
	}

	Deserializer& operator>>(unsigned long long &x)
	{
		read(x);
		x = big_endian_byteswap(x);
		return *this;
	}

	Deserializer& operator>>(double &x)
	{
		read(x);
		return *this;
	}

	Deserializer& operator>>(std::string& s)
	{
		s.clear();
		if (!read_to(std::back_inserter(s), '\0'))
			throw EndOfStream();
		return *this;
	}

	template<typename T>
	void read(T &x)
	{
		if (avail() >= (int64_t)sizeof(T)) {
#ifdef __sparc__
			std::copy(buffer_->begin, buffer_->begin + sizeof(T), (char*)&x);
#else
			x = *(const T*)buffer_->begin;
#endif
			buffer_->begin += sizeof(T);
		}
		else if (read_raw((char*)&x, sizeof(T)) != sizeof(T))
			throw EndOfStream();
	}

	template<class T>
	size_t read(T *ptr, size_t count)
	{
		return read_raw((char*)ptr, count * sizeof(T)) / sizeof(T);
	}

	const char* data() const
	{
		return buffer_->begin;
	}

	template<typename It>
	bool read_to(It dst, char delimiter)
	{
		int d = delimiter;
		do {
			const char* p = (const char*)memchr((void*)buffer_->begin, d, avail());
			if (p == 0) {
				std::copy(buffer_->begin, buffer_->end, dst);
			}
			else {
				const size_t n = p - buffer_->begin;
				std::copy(buffer_->begin, buffer_->begin + n, dst);
				buffer_->begin += n + 1;
				return true;
			}
		} while (buffer_->fetch());
		return false;
	}

	template<typename It>
	std::pair<bool, int64_t> read_to(It dst, char line_delimiter, char record_start)
	{
		int d = line_delimiter;
		int64_t nl = 0;
		do {
			const char* p = (const char*)memchr((void*)buffer_->begin, d, avail());
			if (p == nullptr) {
				std::copy(buffer_->begin, buffer_->end, dst);
				if (!buffer_->fetch())
					return { false,nl };
			}
			else {
				++nl;
				if (p + 1 < buffer_->begin + avail()) {
					if (p[1] == record_start) {
						const size_t n = p - buffer_->begin;
						std::copy(buffer_->begin, buffer_->begin + n, dst);
						buffer_->begin += n + 1;
						return { true,nl };
					}
					else {
						std::copy(buffer_->begin, p + 1, dst);
						buffer_->begin = p + 1;
					}
				}
				else {
					std::copy(buffer_->begin, p, dst);
					if (!buffer_->fetch())
						return { false,nl };
					if (buffer_->begin[0] == record_start)
						return { true,nl };
					else
						*dst++ = line_delimiter;
				}
			}
		} while (true);
	}

	size_t read_raw(char *ptr, int64_t count);
	std::string peek(int64_t n);
	int64_t file_size() {
		return buffer_->file_size();
	}
	~Deserializer();
	FILE* file();

protected:
	
	void pop(char *dst, int64_t n);

	int64_t avail() const
	{
		return buffer_->end - buffer_->begin;
	}
	
	InputStreamBuffer *buffer_;

};
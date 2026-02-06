/****
Copyright © 2012-2026 Benjamin J. Buchfink <buchfink@gmail.com>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

	http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
****/
// SPDX-License-Identifier: Apache-2.0

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
		read(&x);
		x = big_endian_byteswap(x);
		return *this;
	}

	Deserializer& operator>>(int32_t &x)
	{
		read(&x);
		x = big_endian_byteswap(x);
		return *this;
	}

	Deserializer& operator>>(int64_t& x)
	{
		read(&x);
		x = big_endian_byteswap(x);
		return *this;
	}

	Deserializer& operator>>(short& x) {
		read(&x);
		x = big_endian_byteswap(x);
		return *this;
	}

	Deserializer& operator>>(unsigned short& x) {
		read(&x);
		x = big_endian_byteswap(x);
		return *this;
	}

	Deserializer& operator>>(unsigned long &x)
	{
		read(&x);
		x = big_endian_byteswap(x);
		return *this;
	}

	Deserializer& operator>>(unsigned long long &x)
	{
		read(&x);
		x = big_endian_byteswap(x);
		return *this;
	}

	Deserializer& operator>>(double &x)
	{
		read(&x);
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
	void read(T* x)
	{
		if (avail() >= (int64_t)sizeof(T)) {
			std::memcpy((void*)x, buffer_->begin, sizeof(T));
			buffer_->begin += sizeof(T);
		}
		else if (read_raw((char*)x, sizeof(T)) != sizeof(T))
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
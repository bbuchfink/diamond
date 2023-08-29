/****
DIAMOND protein aligner
Copyright (C) 2016-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
						
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
#include <vector>
#include <utility>
#include <iterator>
#include <string.h>
#include <assert.h>
#include "stream_entity.h"
#include "../algo/varint.h"
#include "../system/endianness.h"

struct DynamicRecordReader;

struct Deserializer
{

	enum { VARINT = 1 };

	Deserializer(StreamEntity* buffer);
	void rewind();
	Deserializer& seek(int64_t pos);
	void seek_forward(size_t n);
	bool seek_forward(char delimiter);
	void close();

	Deserializer(const char *begin, const char *end, int flags = 0):
		varint(flags & VARINT),
		buffer_(NULL),
		begin_(begin),
		end_(end)
	{}

	Deserializer& operator>>(unsigned &x)
	{
		if (varint)
			read_varint(*this, x);
		else {
			read(x);
			x = big_endian_byteswap(x);
		}
		return *this;
	}

	Deserializer& operator>>(int32_t &x)
	{
		if (varint) {
			uint32_t i;
			read_varint(*this, i);
			x = (int32_t)i;
		}
		else {
			read(x);
			x = big_endian_byteswap(x);
		}
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

	Deserializer& operator>>(std::string &s)
	{
		s.clear();
		if (!read_to(std::back_inserter(s), '\0'))
			throw EndOfStream();
		return *this;
	}

	Deserializer& operator>>(std::vector<std::string> &v)
	{
		uint32_t n;
		varint = false;
		*this >> n;
		v.clear();
		v.reserve(n);
		std::string s;
		for (uint32_t i = 0; i < n; ++i) {
			*this >> s;
			v.push_back(std::move(s));
		}
		return *this;
	}

	Deserializer& operator>>(std::vector<int32_t> &v)
	{
		int32_t n, x;
		*this >> n;
		v.clear();
		v.reserve(n);
		for (int32_t i = 0; i < n; ++i) {
			*this >> x;
			v.push_back(x);
		}
		return *this;
	}

	template<typename Type1, typename Type2>
	Deserializer& operator>>(std::pair<Type1, Type2>& p) {
		*this >> p.first;
		*this >> p.second;
		return *this;
	}

	template<typename _t>
	void read(_t &x)
	{
		if (avail() >= sizeof(_t)) {
			x = *(const _t*)begin_;
			begin_ += sizeof(_t);
		}
		else if (read_raw((char*)&x, sizeof(_t)) != sizeof(_t))
			throw EndOfStream();
	}

	template<class _t>
	size_t read(_t *ptr, size_t count)
	{
		return read_raw((char*)ptr, count*sizeof(_t)) / sizeof(_t);
	}

	const char* data() const
	{
		return begin_;
	}

	template<typename It>
	bool read_to(It dst, char delimiter)
	{
		int d = delimiter;
		do {
			const char* p = (const char*)memchr((void*)begin_, d, avail());
			if (p == 0) {
				std::copy(begin_, end_, dst);
			}
			else {
				const size_t n = p - begin_;
				std::copy(begin_, begin_ + n, dst);
				begin_ += n + 1;
				return true;
			}
		} while (fetch());
		return false;
	}

	template<typename It>
	void read_raw(size_t count, It out)
	{
		if (count <= avail()) {
			pop(count, out);
			return;
		}
		do {
			const size_t n = std::min(count, avail());
			pop(n, out);
			count -= n;
			if (avail() == 0)
				fetch();
		} while (count > 0 && avail() > 0);
	}

	size_t read_raw(char *ptr, size_t count);
	DynamicRecordReader read_record();
	int64_t file_size() {
		return buffer_->file_size();
	}
	~Deserializer();

	bool varint;

protected:
	
	void pop(char *dst, size_t n);

	template<typename It>
	void pop(size_t n, It out)
	{
		assert(n <= avail());
		std::copy(begin_, begin_ + n, out);
		begin_ += n;
	}

	bool fetch();

	size_t avail() const
	{
		return end_ - begin_;
	}
	
	StreamEntity *buffer_;
	const char *begin_, *end_;
};

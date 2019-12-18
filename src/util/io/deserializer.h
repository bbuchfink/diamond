/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

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

#ifndef DESERIALIZER_H_
#define DESERIALIZER_H_

#include <vector>
#include <utility>
#include "stream_entity.h"
#include "../algo/varint.h"

struct DynamicRecordReader;

struct Deserializer
{

	enum { VARINT = 1 };

	Deserializer(StreamEntity* buffer);
	void rewind();
	Deserializer& seek(size_t pos);
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
		else
			read(x);
		return *this;
	}

	Deserializer& operator>>(int &x)
	{
		read(x);
		return *this;
	}

	Deserializer& operator>>(unsigned long &x)
	{
		read(x);
		return *this;
	}

	Deserializer& operator>>(unsigned long long &x)
	{
		read(x);
		return *this;
	}

	Deserializer& operator>>(double &x)
	{
		read(x);
		return *this;
	}

	Deserializer& operator>>(std::string &s)
	{
		if (!read_until(s, '\0'))
			throw EndOfStream();
		return *this;
	}

	Deserializer& operator>>(std::vector<std::string> &v)
	{
		int n;
		*this >> n;
		v.clear();
		v.reserve(n);
		std::string s;
		for (int i = 0; i < n; ++i) {
			*this >> s;
			v.push_back(std::move(s));
		}
		return *this;
	}

	Deserializer& operator>>(std::vector<unsigned> &v)
	{
		unsigned n, x;
		*this >> n;
		v.clear();
		v.reserve(n);
		for (unsigned i = 0; i < n; ++i) {
			*this >> x;
			v.push_back(x);
		}
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

	size_t read_raw(char *ptr, size_t count);
	bool read_until(std::string &dst, char delimiter);
	bool read_until(std::vector<char> &dst, char delimiter);
	DynamicRecordReader read_record();
	~Deserializer();

	bool varint;

protected:

	template<typename _t>
	bool read_to(_t &container, char delimiter);
	void pop(char *dst, size_t n);
	bool fetch();

	size_t avail() const
	{
		return end_ - begin_;
	}
	
	StreamEntity *buffer_;
	const char *begin_, *end_;
};

#endif
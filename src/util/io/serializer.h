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

#pragma once
#include <string>
#include <set>
#include <vector>
#include "../algo/varint.h"
#include "stream_entity.h"
#include "consumer.h"
#include "../system/endianness.h"

struct Serializer : public Consumer
{

	enum { VARINT = 1 };

	Serializer(StreamEntity *buffer, int flags = 0);

	Serializer& operator<<(int32_t x)
	{
		if (varint_)
			write_varint((uint32_t)x, *this);
		else
			write(big_endian_byteswap(x));
		return *this;
	}

	Serializer& operator<<(long long x)
	{
		write(big_endian_byteswap(x));
		return *this;
	}

	Serializer& operator<<(long x)
	{
		write(big_endian_byteswap(x));
		return *this;
	}

	Serializer& operator<<(unsigned x)
	{
		if (varint_)
			write_varint(x, *this);
		else
			write(big_endian_byteswap(x));
		return *this;
	}

	Serializer& operator<<(unsigned long x)
	{
		write(big_endian_byteswap(x));
		return *this;
	}

	Serializer& operator<<(unsigned long long x)
	{
		write(big_endian_byteswap(x));
		return *this;
	}

	Serializer& operator<<(const double x) {
		write(x);
		return *this;
	}

	Serializer& operator<<(const std::vector<unsigned> &v)
	{
		*this << (unsigned)v.size();
		for (std::vector<unsigned>::const_iterator i = v.begin(); i < v.end(); ++i)
			*this << *i;
		return *this;
	}

	Serializer& operator<<(const std::vector<int32_t> &v)
	{
		*this << (unsigned)v.size();
		for (std::vector<int32_t>::const_iterator i = v.begin(); i < v.end(); ++i)
			*this << *i;
		return *this;
	}

	Serializer& operator<<(const std::set<unsigned> &v)
	{
		*this << (unsigned)v.size();
		for (std::set<unsigned>::const_iterator i = v.begin(); i != v.end(); ++i)
			*this << *i;
		return *this;
	}

	Serializer& operator<<(const std::set<int32_t>& v)
	{
		*this << (unsigned)v.size();
		for (std::set<int32_t>::const_iterator i = v.begin(); i != v.end(); ++i)
			*this << *i;
		return *this;
	}

	Serializer& operator<<(const std::string &s)
	{
		write(s.c_str(), s.length() + 1);
		return *this;
	}

	Serializer& operator<<(const std::vector<std::string> &v)
	{
		varint_ = false;
		*this << (uint32_t)v.size();
		for (std::vector<std::string>::const_iterator i = v.begin(); i < v.end(); ++i)
			*this << *i;
		return *this;
	}

	template<typename Type1, typename Type2>
	Serializer& operator<<(const std::pair<Type1, Type2>& p) {
		*this << p.first;
		*this << p.second;
		return *this;
	}

	template<typename _t>
	void write(const _t &x)
	{
		if (sizeof(_t) <= avail()) {
			*(_t*)next_ = x;
			next_ += sizeof(_t);
		}
		else
			write(&x, 1);
	}

	template<typename _t>
	void write(const _t *ptr, size_t count)
	{
		write_raw((const char*)ptr, count * sizeof(_t));
	}

	template<class _t>
	void write_raw(const std::vector<_t> &v)
	{
		write(v.data(), v.size());
	}

	void set(int flag);
	void unset(int flag);
	void write_raw(const char *ptr, size_t count);
	void seek(int64_t p, int origin = SEEK_SET);
	void rewind();
	size_t tell();
	void close();
	std::string file_name() const;
	FILE* file();
	virtual void consume(const char *ptr, size_t n) override;
	virtual void finalize() override;
	~Serializer();

protected:

	void reset_buffer();
	void flush();

	size_t avail() const
	{
		return end_ - next_;
	}

	StreamEntity *buffer_;
	char *begin_, *next_, *end_;
	bool varint_;

};
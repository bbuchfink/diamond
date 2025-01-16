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
#include "stream_entity.h"
#include "consumer.h"
#include "../system/endianness.h"

struct Serializer : public Consumer
{

	Serializer(StreamEntity *buffer);

	Serializer& operator<<(int32_t x)
	{
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

	Serializer& operator<<(const std::string &s)
	{
		write(s.c_str(), s.length() + 1);
		return *this;
	}

	template<typename T>
	void write(const T &x)
	{
		if (sizeof(T) <= avail()) {
#ifdef __sparc__
			std::copy((char*)&x, (char*)&x + sizeof(T), next_);
#else
			*(T*)next_ = x;
#endif
			next_ += sizeof(T);
		}
		else
			write(&x, 1);
	}

	template<typename T>
	void write(const T *ptr, size_t count)
	{
		write_raw((const char*)ptr, count * sizeof(T));
	}

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

};
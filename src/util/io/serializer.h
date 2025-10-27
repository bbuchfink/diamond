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

	int64_t file_size() {
		return buffer_->file_size();
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
	void flush();

protected:

	void reset_buffer();	

	size_t avail() const
	{
		return end_ - next_;
	}

	StreamEntity *buffer_;
	char *begin_, *next_, *end_;

};
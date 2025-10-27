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
#include <array>
#include <vector>
#include <stdexcept>
#include <stdint.h>
#include <utility>
#include <string.h>
#include <string>
#include "algo/varint.h"

struct BinaryBuffer : public std::vector<char>
{

	struct Iterator
	{
		Iterator(std::vector<char>::const_iterator begin, std::vector<char>::const_iterator end):
			ptr_ (begin),
			end_ (end)
		{ }
		Iterator& operator>>(uint32_t &x)
		{ read(x); return *this; }
		Iterator& operator>>(uint8_t &x)
		{ read(x); return *this; }
		template<typename T>
		void read(T &x)
		{
			check(sizeof(T));
#ifdef __sparc__
			std::copy(ptr_, ptr_ + sizeof(T), (char*)&x);
#else
			x = *(T*)(&*ptr_);
#endif
			ptr_ += sizeof(T);
		}
		template<typename T>
		void read(std::vector<T> &v, size_t count)
		{
			const size_t l = sizeof(T) * count;
			check(l);
			v.resize(count);
			memcpy(v.data(), &*ptr_, l);
			ptr_ += l;
		}
		void read_packed(uint8_t length, uint32_t &dst)
		{
			switch(length) {
			case 0: uint8_t x; read(x); dst = x; break;
			case 1: uint16_t y; read(y); dst = y; break;
			case 2: read(dst);
			}
		}
		void read_packed(uint8_t length, int32_t& dst)
		{
			switch (length) {
			case 0: uint8_t x; read(x); dst = (int32_t)x; break;
			case 1: uint16_t y; read(y); dst = (int32_t)y; break;
			case 2: read(dst);
			}
		}
		void read_varint(uint32_t &dst)
		{
			check(1);
			std::array<char, 5> buf;
			std::copy(ptr_, std::min(ptr_ + 5, end_), buf.begin());
			std::pair<uint32_t, const char*> r = read_varuint32(buf.data());
			const int64_t n = r.second - buf.data();
			check(n);
			ptr_ += n;
			dst = r.first;
		}
		Iterator& operator>>(std::string &dst)
		{
			dst.clear();
			char c;
			while(read(c), c != '\0')
				dst.push_back(c);
			return *this;
		}
		bool good() const
		{ return ptr_ < end_; }
	private:
		void check(size_t size) const
		{
			if (ptr_ + size > end_) throw std::runtime_error("Unexpected end of file.");
		}
		std::vector<char>::const_iterator ptr_, end_;
	};

	Iterator begin() const
	{ return Iterator (std::vector<char>::begin(), std::vector<char>::end()); }

};
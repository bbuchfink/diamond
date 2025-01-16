/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

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
#include <array>
#include <vector>
#include <stdexcept>
#include <stdint.h>
#include <utility>
#include <string.h>
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

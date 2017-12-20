/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#ifndef BINARY_BUFFER_H_
#define BINARY_BUFFER_H_

#include <vector>
#include <stdexcept>
#include <stdint.h>
#include <string.h>
#include "intrin.h"

using std::vector;
using std::string;

struct Binary_buffer : public vector<char>
{

	struct Iterator
	{
		Iterator(vector<char>::const_iterator begin, vector<char>::const_iterator end):
			ptr_ (begin),
			end_ (end)
		{ }
		Iterator& operator>>(uint32_t &x)
		{ read(x); return *this; }
		Iterator& operator>>(uint8_t &x)
		{ read(x); return *this; }
		template<typename _t>
		void read(_t &x)
		{
			check(sizeof(_t));
			x = *(_t*)(&*ptr_);
			ptr_ += sizeof(_t);
		}
		template<typename _t>
		void read(vector<_t> &v, size_t count)
		{
			const size_t l = sizeof(_t) * count;
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
		void read_varint(uint32_t &dst)
		{
			uint8_t b0, b1;
			uint16_t b2;
			uint32_t b3;
			read(b0);
			int c = ctz((uint32_t)b0);
			switch (c) {
			case 0:
				dst = b0 >> 1;
				return;
			case 1:
				read(b1);
				dst = (uint32_t(b1) << 6) | (uint32_t(b0) >> 2);
				return;
			case 2:
				read(b2);
				dst = (uint32_t(b2) << 5) | (uint32_t(b0) >> 3);
				return;
			case 3:
				read(b1);
				read(b2);
				dst = (uint32_t(b2) << 12) | (uint32_t(b1) << 4) | (uint32_t(b0) >> 4);
				return;
			case 4:
				read(b3);
				dst = (b3 << 3) | (uint32_t(b0) >> 5);
				return;
			default:
				throw std::runtime_error("Format error: Invalid varint encoding.");
			}
		}
		Iterator& operator>>(string &dst)
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
		{ if(ptr_+size > end_) throw std::runtime_error("Unexpected end of file."); }
		vector<char>::const_iterator ptr_, end_;
	};

	Iterator begin() const
	{ return Iterator (vector<char>::begin(), vector<char>::end()); }

};

#endif /* BINARY_BUFFER_H_ */

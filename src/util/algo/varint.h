/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

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

#ifndef VARINT_H_
#define VARINT_H_

#include <stdint.h>

template<typename _out>
inline void write_varint(unsigned x, _out &out)
{
	if (x < 1 << 7) {
		out.write((uint8_t)(x << 1 | 1));
	}
	else if (x < 1 << 14) {
		out.write((uint16_t)(x << 2 | 2));
	}
	else if (x < 1 << 21) {
		out.write(uint8_t((x & 31) << 3 | 4));
		out.write(uint16_t(x >> 5));
	}
	else if (x < 1 << 28) {
		out.write(uint32_t(x << 4 | 8));
	}
	else {
		out.write(uint8_t((x & 7) << 5 | 16));
		out.write(uint32_t(x >> 3));
	}
}

#endif
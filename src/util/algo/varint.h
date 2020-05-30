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
#include <stdint.h>
#include <stdexcept>
#include "../intrin.h"
#include "../system/endianness.h"

template<typename _out>
inline void write_varint(uint32_t x, _out &out)
{
	if (x < 1 << 7) {
		out.write((uint8_t)(x << 1 | 1));
	}
	else if (x < 1 << 14) {
		out.write(big_endian_byteswap((uint16_t)(x << 2 | 2)));
	}
	else if (x < 1 << 21) {
		out.write(uint8_t((x & 31) << 3 | 4));
		out.write(big_endian_byteswap(uint16_t(x >> 5)));
	}
	else if (x < 1 << 28) {
		out.write(big_endian_byteswap(uint32_t(x << 4 | 8)));
	}
	else {
		out.write(uint8_t((x & 7) << 5 | 16));
		out.write(big_endian_byteswap(uint32_t(x >> 3)));
	}
}

template<typename _buf>
void read_varint(_buf &buf, uint32_t &dst)
{
	uint8_t b0, b1;
	uint16_t b2;
	uint32_t b3;
	buf.read(b0);
	int c = ctz((uint32_t)b0);
	switch (c) {
	case 0:
		dst = b0 >> 1;
		return;
	case 1:
		buf.read(b1);
		dst = (uint32_t(b1) << 6) | (uint32_t(b0) >> 2);
		return;
	case 2:
		buf.read(b2);
		b2 = big_endian_byteswap(b2);
		dst = (uint32_t(b2) << 5) | (uint32_t(b0) >> 3);
		return;
	case 3:
		buf.read(b1);
		buf.read(b2);
		b2 = big_endian_byteswap(b2);
		dst = (uint32_t(b2) << 12) | (uint32_t(b1) << 4) | (uint32_t(b0) >> 4);
		return;
	case 4:
		buf.read(b3);
		b3 = big_endian_byteswap(b3);
		dst = (b3 << 3) | (uint32_t(b0) >> 5);
		return;
	default:
		throw std::runtime_error("Format error: Invalid varint encoding.");
	}
}
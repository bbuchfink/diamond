/****
DIAMOND protein aligner
Copyright (C) 2016-2024 Max Planck Society for the Advancement of Science e.V.
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

inline char* write_varuint32(uint32_t x, char* out)
{
	if (x < 1 << 7) {
		*(uint8_t*)out = (uint8_t)(x << 1 | 1);
		return out + 1;
	}
	else if (x < 1 << 14) {
		*(uint16_t*)out = big_endian_byteswap((uint16_t)(x << 2 | 2));
		return out + 2;
	}
	else if (x < 1 << 21) {
		*(uint8_t*)out = uint8_t((x & 31) << 3 | 4);
		*(uint16_t*)(out + 1) = big_endian_byteswap(uint16_t(x >> 5));
		return out + 3;
	}
	else if (x < 1 << 28) {
		*(uint32_t*)out = big_endian_byteswap(uint32_t(x << 4 | 8));
		return out + 4;
	}
	else {
		*(uint8_t*)out = uint8_t((x & 7) << 5 | 16);
		*(uint32_t*)(out + 1) = big_endian_byteswap(uint32_t(x >> 3));
		return out + 5;
	}
}

inline std::pair<uint32_t, const char*> read_varuint32(const char* ptr) {
	uint8_t b0, b1;
	uint16_t b2;
	uint32_t b3;
	b0 = *ptr++;
	int c = ctz((uint32_t)b0);
	switch (c) {
	case 0:
		return { b0 >> 1, ptr };
	case 1:
		b1 = *ptr++;
		return { (uint32_t(b1) << 6) | (uint32_t(b0) >> 2),ptr };
	case 2:
		b2 = *(uint16_t*)ptr;
		b2 = big_endian_byteswap(b2);
		return { (uint32_t(b2) << 5) | (uint32_t(b0) >> 3),ptr + 2 };
	case 3:
		b1 = *ptr++;
		b2 = *(uint16_t*)ptr;
		b2 = big_endian_byteswap(b2);
		return { (uint32_t(b2) << 12) | (uint32_t(b1) << 4) | (uint32_t(b0) >> 4),ptr + 2 };
	case 4:
		b3 = *(uint32_t*)ptr;
		b3 = big_endian_byteswap(b3);
		return { (b3 << 3) | (uint32_t(b0) >> 5),ptr + 4 };
	default:
		throw std::runtime_error("Format error: Invalid varint encoding.");
	}
}
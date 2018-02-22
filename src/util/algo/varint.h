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

/*template<typename _t, typename _out>
inline void varint_write_wrapper(const _t &x, _out &out)
{
}*/

template<typename _out>
inline void write_varint(unsigned x, _out &out)
{
	if (x < 1 << 7) {
		varint_write_wrapper((uint8_t)(x << 1 | 1), out);
	}
	else if (x < 1 << 14) {
		varint_write_wrapper((uint16_t)(x << 2 | 2), out);
	}
	else if (x < 1 << 21) {
		varint_write_wrapper(uint8_t((x & 31) << 3 | 4), out);
		varint_write_wrapper(uint16_t(x >> 5), out);
	}
	else if (x < 1 << 28) {
		varint_write_wrapper(uint32_t(x << 4 | 8), out);
	}
	else {
		varint_write_wrapper(uint8_t((x & 7) << 5 | 16), out);
		varint_write_wrapper(uint32_t(x >> 3), out);
	}
}

#endif
/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen
						
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
#include "../util/system.h"

#pragma pack(1)

struct packed_uint40_t
{
	uint8_t		high;
	uint32_t	low;
	packed_uint40_t():
		high (),
		low ()
	{ }
	packed_uint40_t(uint64_t v):
		high ((uint8_t)(v>>32)),
		low ((uint32_t)(v&0xfffffffflu))
	{ }
	packed_uint40_t& operator=(uint32_t x) {
		high = 0;
		low = x;
		return *this;
	}
	operator uint64_t() const
	{ return (uint64_t(high) << 32) | low; }
	operator int64_t() const {
		return (int64_t(high) << 32) | (int64_t)low;
	}
	operator uint32_t() const {
		return low;
	}
	operator int32_t() const {
		return low;
	}
	bool operator==(const packed_uint40_t& rhs) const {
		return high == rhs.high && low == rhs.low;
	}
	bool operator!=(const packed_uint40_t& rhs) const {
		return high != rhs.high && low != rhs.low;
	}
	bool operator<(const packed_uint40_t &rhs) const
	{ return high < rhs.high || (high == rhs.high && low < rhs.low); }
	friend uint64_t operator-(const packed_uint40_t &x, const packed_uint40_t &y)
	{ return (uint64_t)(x) - (uint64_t)(y); }
} PACKED_ATTRIBUTE ;

typedef packed_uint40_t PackedLoc;

#pragma pack()


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

#ifndef PACKED_LOC_H_
#define PACKED_LOC_H_

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
	operator const uint64_t() const
	{ return (uint64_t(high) << 32) | low; }
	bool operator<(const packed_uint40_t &rhs) const
	{ return high < rhs.high || (high == rhs.high && low < rhs.low); }
	friend uint64_t operator-(const packed_uint40_t &x, const packed_uint40_t &y)
	{ return (const uint64_t)(x) - (const uint64_t)(y); }
} PACKED_ATTRIBUTE ;

typedef packed_uint40_t Packed_loc;
typedef size_t Loc;

#pragma pack()

#endif /* PACKED_LOC_H_ */

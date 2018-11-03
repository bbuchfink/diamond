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

#ifndef INTRIN_H_
#define INTRIN_H_

#include <stdint.h>
#ifdef _MSC_VER
#include <intrin.h>
#endif

inline unsigned popcount32(unsigned x)
{
#ifdef _MSC_VER
	return __popcnt(x);
#else
	return __builtin_popcount(x);
#endif
}

inline unsigned popcount64(unsigned long long x)
{
#ifdef _MSC_VER
	return (unsigned)__popcnt64(x);
#else
	return __builtin_popcountll(x);
#endif
}

inline int ctz(uint32_t x)
{
#ifdef _MSC_VER
	unsigned long i;
	unsigned char c = _BitScanForward(&i, x);
	return i;
#else
	return __builtin_ctz(x);
#endif
}

inline int ctz(uint64_t x)
{
#ifdef _MSC_VER
	if (x)
		return (int)__popcnt64((x ^ (x - 1)) >> 1);
	else
		return CHAR_BIT * sizeof(x);
#else
	return __builtin_ctzll(x);
#endif
}

#endif
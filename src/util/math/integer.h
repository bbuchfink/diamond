/****
DIAMOND protein aligner
Copyright (C) 2020 Max Planck Society for the Advancement of Science e.V.

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
#include <limits.h>
#include <stdint.h>
#include <algorithm>
#include <math.h>

static inline int8_t saturated_add(int8_t x, int8_t y) {
	return (int8_t)std::max(int32_t(x) + int32_t(y), (int32_t)std::numeric_limits<int8_t>::min());
}

static inline int16_t saturated_add(int16_t x, int16_t y) {
	return (int16_t)std::max(int32_t(x) + int32_t(y), (int32_t)std::numeric_limits<int16_t>::min());
}

static inline int32_t saturated_add(int32_t x, int32_t y) {
	return x + y;
}

static inline size_t bit_length(size_t x) {
	return (size_t)ceil(log(x) / log(2)) + 1;
}

static inline uint64_t next_power_of_2(double x)
{
	return 1llu << uint64_t(ceil(log(x) / log(2)));
}

template<typename _it, int N>
bool next_combination(_it begin, _it end) {
	_it i = begin;
	while (i < end) {
		if (*i < N - 1) {
			++(*i);
			return true;
		}
		else {
			*i = 0;
			++i;
		}
	}
	return false;
}

template<typename I>
I power(I x, I p)
{
	if (p == 0) return 1;
	if (p == 1) return x;

	const I t = power(x, p / 2);
	if (p % 2 == 0)
		return t * t;
	else
		return x * t * t;
}
/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

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

#ifndef RADIX_SORT_H_
#define RADIX_SORT_H_

#include <utility>
#include <stdint.h>
#include "radix_cluster.h"
#include "../math/integer.h"

template<typename _t>
void radix_sort(_t* begin, _t* end, uint32_t max_key) {
	const size_t n = end - begin;
	if (n <= 1)
		return;
	const size_t bit_len = bit_length(max_key), rounds = (bit_len + config.radix_bits - 1) / config.radix_bits;
	_t* buf = (_t*)new char[n * sizeof(_t)];

	_t* in = begin, * out = buf;
	for (int i = 0; i < rounds; ++i) {
		parallel_radix_cluster(Relation<_t>(in, n), i * config.radix_bits, out);
		std::swap(in, out);
	}

	if (out == begin)
		std::copy(buf, buf + n, begin);

	delete[] buf;
}

#endif
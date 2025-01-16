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
#include <stdint.h>
#include "radix_cluster.h"
#include "../math/integer.h"

template<typename T, typename GetKey>
void radix_sort(T* begin, T* end, uint32_t max_key, size_t threads) {
	typedef typename T::Key Key;
	const size_t n = end - begin;
	if (n <= 1)
		return;
	const uint32_t bit_len = (uint32_t)bit_length(max_key), rounds = (bit_len + config.radix_bits - 1) / config.radix_bits;
	T* buf = (T*)new char[n * sizeof(T)];

	T* in = begin, * out = buf;
	for (uint32_t i = 0; i < rounds; ++i) {
		if(threads > 1)
			parallel_radix_cluster<T, GetKey>(Relation<T>(in, n), i * config.radix_bits, out, threads);
		else {
			unsigned *hst = new unsigned[(size_t)1 << config.radix_bits];
			radix_cluster<T, GetKey>(Relation<T>(in, n), i * config.radix_bits, out, hst);
			delete[] hst;
		}

		std::swap(in, out);
	}

	if (out == begin)
		std::copy(buf, buf + n, begin);

	delete[] buf;
}
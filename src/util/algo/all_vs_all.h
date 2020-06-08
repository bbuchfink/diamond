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
#include <algorithm>
#include "../data_structures/flat_array.h"

template<typename _t>
void all_vs_all(const _t* a, uint32_t na, const _t* b, uint32_t nb, FlatArray<uint32_t>& out) {
	for (uint32_t i = 0; i < na; ++i) {
		const _t e = a[i];
		out.next();
		for (uint32_t j = 0; j < nb; ++j)
			if (e == b[j])
				out.push_back(j);
	}
}

template<typename _t, typename _f>
void all_vs_all(const _t* a, uint32_t na, const _t* b, uint32_t nb, uint32_t tile_size, _f &callback) {
	thread_local FlatArray<uint32_t> out;
	for (uint32_t i = 0; i < na; i += tile_size) {
		for (uint32_t j = 0; j < nb; j += tile_size) {
			out.clear();
			all_vs_all(a + i, std::min(tile_size, na - i), b + j, std::min(tile_size, nb - j), out);
			callback(out, i, j);
		}
	}
}
/****
DIAMOND protein aligner
Copyright (C) 2022 Max Planck Society for the Advancement of Science e.V.

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
#include <utility>
#include <iterator>
#include "diagonal_segment.h"

namespace Util { namespace Geometry {

using Hit = std::pair<Loc, Loc>;

inline bool cmp_diag(const Hit& a, const Hit& b) {
	const Loc d1 = a.first - a.second, d2 = b.first - b.second;
	return d1 < d2 || (d1 == d2 && a.first < b.first);
}

// Merges seed hits on one diagonal assumed to be sorted in ascending order.
template<typename It>
inline std::vector<DiagonalSegment> merge_hits(const It begin, const It end, int32_t kmer_size, Loc window, Loc min_len) {
	std::vector<DiagonalSegment> v;
	if (begin == end)
		return v;
	Loc a = begin->first, b = begin->first + kmer_size, d = begin->first - begin->second;
	for (It i = std::next(begin); i != end; ++i) {
		if (i->first - b < window)
			b = std::max(b, i->first + kmer_size);
		else {
			if (b - a >= min_len)
				v.emplace_back(a, a - d, b - a, 0);
			a = i->first;
			b = i->first + kmer_size;
		}
	}
	if (b - a >= min_len)
		v.emplace_back(a, a - d, b - a, 0);
	return v;
}

}}

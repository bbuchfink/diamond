/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once
#include <vector>
#include "table.h"

using std::vector;

namespace Util { namespace Tsv {

struct BuildHelper {
	BuildHelper(int64_t size) {
		static const double TOLERANCE = 1.1;
		buffer.reserve(size_t(size * TOLERANCE));
		limits.push_back(0);
	}
	void add(Table* t) {
		buffer.insert(buffer.end(), t->data_.begin(), t->data_.end());
		limits.insert(limits.end(), t->limits_.begin() + 1, t->limits_.end());
		counts.push_back((int32_t)t->size());
	}
	Table get(const Schema& schema) {
		const int64_t n = counts.size();
		if (n == 0)
			return Table(schema);
		vector<int64_t> offsets;
		offsets.reserve(n);
		offsets.push_back(0);
		vector<int64_t>::const_iterator it = limits.begin() + counts[0];
		for (int64_t i = 1; i < n; ++i) {
			offsets.push_back(*it + offsets.back());
			it += counts[i];
		}
		vector<int64_t>::iterator begin = limits.begin() + counts[0] + 1;
		for (int64_t i = 1; i < n; ++i) {
			vector<int64_t>::iterator end = begin + counts[i];
			const auto d = offsets[i];
			for (auto it = begin; it < end; ++it)
				*it += d;
			begin += counts[i];
		}
		return Table(schema, std::move(buffer), std::move(limits));
	}
	vector<char> buffer;
	vector<int64_t> limits;
	vector<int32_t> counts;
};

}}
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

#include <queue>
#include "file.h"

using std::priority_queue;
using std::pair;
using std::vector;
using std::greater;

namespace Util { namespace Tsv {

File* merge(std::vector<File*>::iterator begin, std::vector<File*>::iterator end, int column) {
	using Pair = pair<int64_t, int64_t>;
	File* out = new File((*begin)->schema(), "", Flags::TEMP);
	priority_queue<Pair, vector<Pair>, greater<Pair>> queue;
	vector<Table> tables;
	tables.reserve(end - begin);
	for (auto it = begin; it != end; ++it) {
		tables.push_back((*it)->read_record());
		if (tables.back().empty())
			continue;
		queue.emplace(tables.back().front().get<int64_t>(column), it - begin);
	}
	while (!queue.empty()) {
		const int64_t f = queue.top().second;
		out->write(tables[f].front());
		queue.pop();
		tables[f] = begin[f]->read_record();
		if (!tables[f].empty())
			queue.emplace(tables[f].front().get<int64_t>(column), f);
	}
	return out;
}

}}
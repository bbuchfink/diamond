/****
Copyright (C) 2012-2026 Benjamin J. Buchfink <buchfink@gmail.com>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

	http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
****/
// SPDX-License-Identifier: Apache-2.0

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
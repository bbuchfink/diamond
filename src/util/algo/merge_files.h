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

template<typename T, typename It, typename F>
void merge_sorted_files(const It begin, const It end, F& f) {
	struct Entry {
		Entry(T&& value, ptrdiff_t idx) :
			value(value),
			idx(idx)
		{
		}
		Entry(const Entry& e) :
			value(e.value),
			idx(e.idx)
		{
		}
		Entry& operator=(const Entry& e) {
			value = e.value;
			idx = e.idx;
			return *this;
		}
		bool operator<(const Entry& e) const {
			return !(value < e.value);
		}
		T value;
		ptrdiff_t idx;
	};
	std::priority_queue<Entry> q;
	std::vector<File*> d;
	T h;
	for (It i = begin; i != end; ++i) {
		d.emplace_back(*i);
		if (deserialize(d.back(), h))
			q.emplace(std::move(h), i - begin);
	}
	while (!q.empty()) {
		f(q.top().value);
		const ptrdiff_t idx = q.top().idx;
		q.pop();
		if (deserialize(d[idx], h))
			q.emplace(std::move(h), idx);
	}
}
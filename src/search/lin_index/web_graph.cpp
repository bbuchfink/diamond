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

#include <algorithm>
#include <stdexcept>
#include "web_graph.h"

using std::vector;

namespace Search {

WebGraph::WebGraph(const Vertex vertex_count) :
	offsets_((size_t)vertex_count + 1, 0)
{}

bool WebGraph::contains(Vertex a, Vertex b) const {
	if (a < b)
		std::swap(a, b);
	if (a >= vertex_count())
		return false;
	const auto begin = neighbors_.begin() + offsets_[a], end = neighbors_.begin() + offsets_[a + 1];
	return std::binary_search(begin, end, b);
}

void WebGraph::insert(vector<Edge>& edges) {
	const Vertex n = vertex_count();
	// Orient every edge towards its larger endpoint, whose list it is stored in.
	for (Edge& e : edges) {
		if (e.first < e.second)
			std::swap(e.first, e.second);
		if (e.first >= n)
			throw std::out_of_range("WebGraph::insert");
	}
	std::sort(edges.begin(), edges.end());
	edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
	if (edges.empty())
		return;

	vector<int64_t> offsets((size_t)n + 1);
	vector<Vertex> neighbors;
	neighbors.reserve(neighbors_.size() + edges.size());
	auto e = edges.cbegin();
	for (Vertex u = 0; u < n; ++u) {
		offsets[u] = (int64_t)neighbors.size();
		auto old = neighbors_.cbegin() + offsets_[u];
		const auto old_end = neighbors_.cbegin() + offsets_[u + 1];
		// Merge the sorted old list of u with its sorted new edges.
		while (e != edges.cend() && e->first == u) {
			while (old != old_end && *old < e->second)
				neighbors.push_back(*old++);
			if (old != old_end && *old == e->second)
				++old;
			neighbors.push_back(e->second);
			++e;
		}
		neighbors.insert(neighbors.end(), old, old_end);
	}
	offsets[n] = (int64_t)neighbors.size();
	offsets_ = std::move(offsets);
	neighbors_ = std::move(neighbors);
}

}

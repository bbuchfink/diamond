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
#include <stdint.h>
#include <utility>
#include <vector>

namespace Search {

/* Undirected graph on a fixed set of vertices numbered from 0, held in memory as
   adjacency lists in compressed sparse row form. Every edge is stored once, in the list
   of its larger endpoint, whose neighbors are kept sorted, so that testing an edge is a
   binary search in a single list.

   The graph is built up in batches: insert merges a batch of edges into the lists and
   rebuilds them, and is therefore meant to be called rarely and with many edges at a
   time. Lookups do not modify the graph and may run concurrently, but not concurrently
   with an insert. */
struct WebGraph {

	using Vertex = uint32_t;
	using Edge = std::pair<Vertex, Vertex>;

	explicit WebGraph(const Vertex vertex_count);

	bool contains(Vertex a, Vertex b) const;
	/* Merges the edges into the graph. Their direction does not matter, and edges that
	   occur more than once or are already part of the graph are only stored once. The
	   vector is reordered. */
	void insert(std::vector<Edge>& edges);

	Vertex vertex_count() const {
		return (Vertex)(offsets_.size() - 1);
	}

	int64_t edge_count() const {
		return (int64_t)neighbors_.size();
	}

	// Bytes of the adjacency lists, i.e. of the offsets and the neighbors stored.
	int64_t data_size() const {
		return (int64_t)(offsets_.size() * sizeof(int64_t) + neighbors_.size() * sizeof(Vertex));
	}

private:

	// Begin of the adjacency list of each vertex in neighbors_, followed by the end of the
	// last one.
	std::vector<int64_t> offsets_;
	std::vector<Vertex> neighbors_;

};

}

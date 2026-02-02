/****
Copyright © 2013-2026 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-2-Clause

#include <queue>
#include "algo.h"
#include "../log_stream.h"

using std::vector;
using std::priority_queue;
using std::pair;
using std::numeric_limits;
using std::swap;
using std::queue;

namespace Util { namespace Algo {

template<typename Int, typename It>
static Int neighbor_count2(It begin, It end, const vector<Int>& centroids) {
	static constexpr Int NIL = std::numeric_limits<Int>::max();
	Int n = 0, last = NIL;
	for (It i = begin; i != end; ++i)
		if (centroids[i->node2] == NIL && i->node2 != last) {
			++n;
			last = i->node2;
		}
	return n;
}

template<typename Int, typename It>
static Int neighbor_count(It begin, It end, const vector<Int>& centroids) {
	static constexpr Int NIL = std::numeric_limits<Int>::max();
	Int n = 0, last = NIL;
	It w = begin;
	for (It i = begin; i != end; ++i) {
		if (i->node1 == NIL)
			break;
		if (centroids[i->node2] == NIL && i->node2 != last) {
			++n;
			last = i->node2;
			if (w < i)
				swap(w, i);
			++w;
		}
	}
	if (w < end)
		w->node1 = NIL;
	return n;
}

template<typename Int, typename It>
static Int neighbor_count(Int node, It begin, It end, const vector<Int>& centroids, const Int* member_counts) {
	static constexpr Int NIL = std::numeric_limits<Int>::max();
	Int n = member_counts[node];
	for (It i = begin; i != end; ++i)
		if (centroids[i->node2] == NIL)
			n += member_counts[i->node2];
	return n;
}

template<typename Int>
static void fix_assignment(vector<Int>& centroids) {
	for (Int i = 0; i < (Int)centroids.size();) {
		if (centroids[centroids[i]] != centroids[i])
			centroids[i] = centroids[centroids[i]];
		else
			++i;
	}
}

template<typename Int>
void make_cluster_gvc(Int rep, FlatArray<Edge<Int>>& neighbors, vector<Int>& centroids, bool merge_recursive) {
	static constexpr Int NIL = std::numeric_limits<Int>::max();
	centroids[rep] = rep;
	for (auto i = neighbors.cbegin(rep); i != neighbors.cend(rep); ++i)
		if (centroids[i->node2] == NIL || (merge_recursive && centroids[i->node2] == i->node2))
			centroids[i->node2] = rep;
}

template<typename Int>
void make_cluster_cc(Int rep, FlatArray<Edge<Int>>& neighbors, vector<Int>& centroids, Int depth) {
	static constexpr Int NIL = std::numeric_limits<Int>::max();
	struct Entry {
		Entry(Int node, Int depth):
			node(node),
			depth(depth)
		{}
		Int node;
		Int depth;
	};
	centroids[rep] = rep;
	queue<Entry> q;
	for (auto i = neighbors.cbegin(rep); i != neighbors.cend(rep); ++i)
		if (centroids[i->node2] == NIL)
			q.emplace(i->node2, 1);
	while (!q.empty()) {
		const Entry node = q.front();
		q.pop();
		if (centroids[node.node] != NIL || node.depth > depth)
			continue;
		for (auto i = neighbors.cbegin(node.node); i != neighbors.cend(node.node); ++i)
			if (centroids[i->node2] == NIL)
				q.emplace(i->node2, node.depth + 1);
		centroids[node.node] = rep;
	}
}

template<typename Int>
vector<Int> greedy_vertex_cover(FlatArray<Edge<Int>>& neighbors, const Int* member_counts, bool merge_recursive, bool reassign, Int connected_component_depth) {
	static constexpr Int NIL = std::numeric_limits<Int>::max();
	TaskTimer timer("Computing edge counts");
	priority_queue<pair<Int, Int>> q;
	vector<Int> centroids(neighbors.size(), NIL);
	for (Int i = 0; i < neighbors.size(); ++i)
		q.emplace(member_counts ? neighbor_count(i, neighbors.cbegin(i), neighbors.cend(i), centroids, member_counts) :
			(Int)neighbors.count(i), i);

	timer.go("Computing vertex cover");
	int64_t cluster_count = 0;
	while (!q.empty()) {
		const Int node = q.top().second;
		q.pop();
		if (centroids[node] != NIL)
			continue;
		const Int count = member_counts ? neighbor_count(node, neighbors.cbegin(node), neighbors.cend(node), centroids, member_counts) :
			neighbor_count(neighbors.begin(node), neighbors.end(node), centroids);
		if (!q.empty() && count < q.top().first)
			q.emplace(count, node);
		else {
			if (connected_component_depth > 0)
				make_cluster_cc(node, neighbors, centroids, connected_component_depth);
			else
				make_cluster_gvc(node, neighbors, centroids, merge_recursive);
			++cluster_count;
		}
	}
	timer.finish();
	message_stream << "Cluster count = " << cluster_count << std::endl;

	if (reassign) {
		timer.go("Computing reassignment");
		vector<double> weights(neighbors.size(), numeric_limits<double>::lowest());
		for (Int node = 0; node < neighbors.size(); ++node)
			if (centroids[node] == node)
				for (auto i = neighbors.cbegin(node); i != neighbors.cend(node); ++i)
					if (centroids[i->node2] != i->node2 && i->weight > weights[i->node2]) {
						weights[i->node2] = i->weight;
						centroids[i->node2] = node;
					}
	}

	if (merge_recursive) {
		timer.go("Computing merges");
		fix_assignment(centroids);
	}

	return centroids;
}

template vector<uint32_t> greedy_vertex_cover<uint32_t>(FlatArray<Edge<uint32_t>>&, const uint32_t*, bool, bool, uint32_t);
template vector<uint64_t> greedy_vertex_cover<uint64_t>(FlatArray<Edge<uint64_t>>&, const uint64_t*, bool, bool, uint64_t);

}}
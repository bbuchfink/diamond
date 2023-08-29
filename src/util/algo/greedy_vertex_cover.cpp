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

#include <algorithm>
#include <float.h>
#include <queue>
#include "algo.h"
#include "../log_stream.h"

using std::vector;
using std::priority_queue;
using std::pair;
using std::numeric_limits;
using std::swap;

namespace Util { namespace Algo {

template<typename Int, typename It>
static Int neighbor_count2(It begin, It end, const vector<Int>& centroids) {
	Int n = 0, last = -1;
	for (It i = begin; i != end; ++i)
		if (centroids[i->node2] == -1 && i->node2 != last) {
			++n;
			last = i->node2;
		}
	return n;
}

template<typename Int, typename It>
static Int neighbor_count(It begin, It end, const vector<Int>& centroids) {
	Int n = 0, last = -1;
	It w = begin;
	for (It i = begin; i != end; ++i) {
		if (i->node1 == -1)
			break;
		if (centroids[i->node2] == -1 && i->node2 != last) {
			++n;
			last = i->node2;
			if (w < i)
				swap(w, i);
			++w;
		}
	}
	if (w < end)
		w->node1 = -1;
	return n;
}

template<typename Int, typename It>
static Int neighbor_count(Int node, It begin, It end, const vector<Int>& centroids, const SuperBlockId* member_counts) {
	Int n = member_counts[node];
	for (It i = begin; i != end; ++i)
		if (centroids[i->node2] == -1)
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
vector<Int> greedy_vertex_cover(FlatArray<Edge<Int>>& neighbors, const SuperBlockId* member_counts, bool merge_recursive) {
	TaskTimer timer("Computing edge counts");
	priority_queue<pair<Int, Int>> q;
	vector<Int> centroids(neighbors.size(), -1);
	for (Int i = 0; i < neighbors.size(); ++i)
		q.emplace(member_counts ? neighbor_count(i, neighbors.cbegin(i), neighbors.cend(i), centroids, member_counts) :
			(Int)neighbors.count(i), i);

	timer.go("Computing vertex cover");
	while (!q.empty()) {
		const Int node = q.top().second;
		q.pop();
		if (centroids[node] != -1)
			continue;
		const Int count = member_counts ? neighbor_count(node, neighbors.cbegin(node), neighbors.cend(node), centroids, member_counts) :
			neighbor_count(neighbors.begin(node), neighbors.end(node), centroids);
		if (!q.empty() && count < q.top().first)
			q.emplace(count, node);
		else {
			centroids[node] = node;
			for (auto i = neighbors.cbegin(node); i != neighbors.cend(node); ++i)
				if (centroids[i->node2] == -1 || (merge_recursive && centroids[i->node2] == i->node2))
					centroids[i->node2] = node;
		}
	}

	timer.go("Computing reassignment");
	vector<double> weights(neighbors.size(), numeric_limits<double>::lowest());
	for (Int node = 0; node < neighbors.size(); ++node)
		if (centroids[node] == node)
			for (auto i = neighbors.cbegin(node); i != neighbors.cend(node); ++i)
				if (centroids[i->node2] != i->node2 && i->weight > weights[i->node2]) {
					weights[i->node2] = i->weight;
					centroids[i->node2] = node;
				}

	if (merge_recursive) {
		timer.go("Computing merges");
		fix_assignment(centroids);
	}

	return centroids;
}

template vector<int32_t> greedy_vertex_cover<int32_t>(FlatArray<Edge<int32_t>>& neighbors, const SuperBlockId* member_counts, bool merge_recursive);
template vector<int64_t> greedy_vertex_cover<int64_t>(FlatArray<Edge<int64_t>>& neighbors, const SuperBlockId* member_counts, bool merge_recursive);

}}
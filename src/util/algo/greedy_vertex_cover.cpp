#include <algorithm>
#include <iterator>
#include <map>
#include "algo.h"

using namespace std;

namespace Util { namespace Algo {

struct GreedyVertexCover {

	GreedyVertexCover(vector<vector<int>> &neighbors) :
		reverse_neighbors(neighbors.size()) {
		const int n = (int)neighbors.size();
		centroid.insert(centroid.begin(), n, -1);
		for (int i = 0; i < n; ++i) {
			for (int j : neighbors[i])
				reverse_neighbors[j].push_back(i);
		}

		vector<int> merged;
		for (int i = 0; i < n; ++i) {
			std::sort(neighbors[i].begin(), neighbors[i].end());
			std::sort(reverse_neighbors[i].begin(), reverse_neighbors[i].end());
			merged.clear();
			std::set_union(neighbors[i].begin(), neighbors[i].end(), reverse_neighbors[i].begin(), reverse_neighbors[i].end(), std::back_inserter(merged));
			neighbors[i] = merged;
		}

		its.reserve(n);
		for (int i = 0; i < n; ++i) {
			its.push_back(count_to_idx.insert(make_pair((int)neighbors[i].size(), i)));
		}

		while (!count_to_idx.empty()) {
			const int i = count_to_idx.rbegin()->second;
			assign_centroid(i, i, neighbors);
			for (int j : neighbors[i])
				if (centroid[j] == -1)
					assign_centroid(j, i, neighbors);
		}
	}

	void assign_centroid(int i, int c, const vector<vector<int>>& neighbors) {
		centroid[i] = c;
		count_to_idx.erase(its[i]);
		for (int j : neighbors[i]) {
			if (centroid[j] >= 0)
				continue;
			auto it = count_to_idx.insert(make_pair(its[j]->first - 1, j));
			count_to_idx.erase(its[j]);
			its[j] = it;
		}
	}

	vector<int> centroid;
	multimap<int, int> count_to_idx;
	vector<multimap<int, int>::iterator> its;
	vector<vector<int>> reverse_neighbors;

};

vector<int> greedy_vertex_cover(vector<vector<int>> &neighbors) {
	return GreedyVertexCover(neighbors).centroid;
}

}}
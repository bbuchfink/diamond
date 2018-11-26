#include <algorithm>
#include <map>
#include "algo.h"

using namespace std;

struct GreedyVortexCover {

	GreedyVortexCover(const vector<vector<int>> &neighbors) :
		reverse_neighbors(neighbors.size()) {
		const int n = (int)neighbors.size();
		centroid.insert(centroid.begin(), n, -1);
		its.reserve(n);
		for (int i = 0; i < n; ++i) {
			its.push_back(count_to_idx.insert(make_pair((int)neighbors[i].size(), i)));
			for (int j : neighbors[i]) {
				reverse_neighbors[j].push_back(i);
			}
		}

		while (!count_to_idx.empty()) {
			const int i = count_to_idx.rbegin()->second;
			assign_centroid(i, i);
			for (int j : neighbors[i])
				if (centroid[j] == -1)
					assign_centroid(j, i);
		}
	}

	void assign_centroid(int i, int c) {
		centroid[i] = c;
		count_to_idx.erase(its[i]);
		for (int j : reverse_neighbors[i]) {
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

vector<int> greedy_vortex_cover(const vector<vector<int>> &neighbors) {
	return GreedyVortexCover(neighbors).centroid;
}

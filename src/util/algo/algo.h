#ifndef UTIL_ALGO_ALGO_H_
#define UTIL_ALGO_ALGO_H_

#include <vector>

namespace Util { namespace Algo {

struct Edge {
	int v1, v2, weight;
	bool operator<(const Edge &x) const {
		return weight > x.weight;
	}
};

std::vector<int> greedy_vortex_cover(const std::vector<std::vector<int>> &neighbors);
std::vector<int> greedy_vortex_cover_weighted(std::vector<Edge> &edges, int vortex_count);

}}

#endif
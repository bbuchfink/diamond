#ifndef UTIL_ALGO_ALGO_H_
#define UTIL_ALGO_ALGO_H_

#include <stdint.h>
#include <vector>

namespace Util { namespace Algo {

struct Edge {
	uint32_t v1, v2, weight;
	bool operator<(const Edge &x) const {
		return weight > x.weight;
	}
};

std::vector<int> greedy_vortex_cover(const std::vector<std::vector<int>> &neighbors);

}}

#endif
#pragma once
#include <stdint.h>
#include <vector>
#include "partition.h"

namespace Util { namespace Algo {

struct Edge {
	uint32_t v1, v2, weight;
	bool operator<(const Edge &x) const {
		return weight > x.weight;
	}
};

std::vector<int> greedy_vertex_cover(std::vector<std::vector<int>> &neighbors);

template<typename It, typename Out>
size_t merge_capped(It i0, const It i1, It j0, const It j1, const size_t cap, Out out) {
	const ptrdiff_t m = (ptrdiff_t)cap;
	ptrdiff_t n = 0;
	size_t count = 0;
	while (n < m) {
		if (i0 == i1) {
			const ptrdiff_t d = std::min(m - n, j1 - j0);
			std::copy(j0, j0 + d, out);
			return count + (size_t)d;
		}
		if (j0 == j1) {
			std::copy(i0, i0 + std::min(m - n, i1 - i0), out);
			return count;
		}
		if (*i0 < *j0)
			*out++ = *i0++;
		else {
			*out++ = *j0++;
			++count;
		}
		++n;
	}
	return count;
}

template<typename It, typename Key>
std::vector<It> partition_table(It begin, It end, size_t n, Key key) {
	std::vector<It> v;
	const size_t count = size_t(end - begin);
	if (count == 0)
		return v;
	Partition<size_t> p(count, n);
	v.reserve(p.parts);
	It e = begin;
	v.push_back(e);
	for (size_t i = 0; i < p.parts; ++i) {
		It it = begin + p.end(i);
		if (it <= e)
			continue;
		const auto k = key(*(it - 1));
		while (it < end && key(*it) == k) ++it;
		v.push_back(it);
		e = it;
	}
	return v;
}

}}
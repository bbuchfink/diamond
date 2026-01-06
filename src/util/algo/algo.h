#pragma once
#include <stdint.h>
#include <numeric>
#include <algorithm>
#include <vector>
#include "partition.h"
#define _REENTRANT
#include "../data_structures/flat_array.h"
#include "basic/value.h"

namespace Util { namespace Algo {

template<typename Int>
struct Edge {
	Edge()
	{}
	using Key = Int;
	Edge(Key node1, Key node2, double weight) :
		node1(node1),
		node2(node2),
		weight(weight)
	{}
	bool operator<(const Edge& e) const {
		return node1 < e.node1 || (node1 == e.node1 && node2 < e.node2);
	}
	struct GetKey {
		Key operator()(const Edge& e) const {
			return e.node1;
		}
	};
	Key node1, node2;
	double weight;
};

template<typename Int>
std::vector<Int> greedy_vertex_cover(FlatArray<Edge<Int>>& neighbors, const Int* member_counts = nullptr, bool merge_recursive = false, bool reassign = true, Int connected_component_depth = 0);
template<typename Int>
std::vector<Int> cluster_pr(FlatArray<Edge<Int>>& neighbors);
//template<typename Int> FlatArray<Int> mcl(FlatArray<Int>& neighbors);

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

template<typename It>
std::vector<std::pair<typename It::value_type, int64_t>> sort_by_value(const It begin, const It end, const int threads) {
	std::vector<std::pair<typename It::value_type, int64_t>> out;
	out.reserve(end - begin);
	for (auto i = begin; i != end; ++i)
		out.emplace_back(*i, i - begin);
	ips4o::parallel::sort(out.begin(), out.end(), std::less<std::pair<typename It::value_type, int64_t>>(), threads);
	return out;
}

struct DSU {
	std::vector<int> parent, rank;
	DSU(int n) {
		parent.resize(n);
		rank.resize(n, 0);
		std::iota(parent.begin(), parent.end(), 0);
	}
	int find(int x) {
		while (x != parent[x]) {
			parent[x] = parent[parent[x]];
			x = parent[x];
		}
		return x;
	}
	void unite(int a, int b) {
		a = find(a);
		b = find(b);
		if (a != b) {
			if (rank[a] < rank[b]) std::swap(a, b);
			parent[b] = a;
			if (rank[a] == rank[b]) rank[a]++;
		}
	}
};


}}
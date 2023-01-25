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

#pragma once
#include <vector>
#include "flat_array.h"

template<typename Int>
struct DisjointSet {

	DisjointSet(Int size) {
		nodes_.reserve(size);
		for (Int i = 0; i < size; ++i)
			nodes_.emplace_back(i, 1);
	}

	Int find(Int i) {
		if (nodes_[i].parent != i)
			return nodes_[i].parent = find(nodes_[i].parent);
		else
			return i;
	}

	void merge(Int i, Int j) {
		i = find(i);
		j = find(j);
		if (i == j)
			return;
		if (nodes_[i].size < nodes_[j].size)
			std::swap(i, j);
		nodes_[j].parent = i;
		nodes_[i].size += nodes_[j].size;
	}

	FlatArray<Int> sets(int threads) {
		std::vector<std::pair<Int, Int>> v;
		v.reserve(nodes_.size());
		for (Int i = 0; i < nodes_.size(); ++i)
			v.emplace_back(find(i), i);
		return make_flat_array(v.begin(), v.end(), threads).first;
	}

private:

	struct Node {
		Node(Int parent, Int size) :
			parent(parent), size(size)
		{}
		Int parent, size;
	};

	std::vector<Node> nodes_;

};
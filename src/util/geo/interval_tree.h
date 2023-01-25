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
#include "../../lib/interval_tree/interval_tree.hpp"
#include "interval.h"

template<typename Int>
struct IntervalTree {

	bool is_overlapped(Interval i, double overlap) {
		Overlaps o(i, Int(i.length() * overlap));
		tree_.overlap_find_all({ i.begin_, i.end_ }, o);
		return o.n >= o.max;
	}

	void insert(Interval i) {
		tree_.insert_overlap({ i.begin_, i.end_ });
	}

private:

	struct Overlaps {
		Overlaps(Interval i, Int max) :
			i(i),
			max(max),
			n(0)
		{}
		bool operator()(typename lib_interval_tree::interval_tree_t<Int>::iterator it) {
			n += intersect(Interval(it->interval().low(), it->interval().high()), i).length();
			return n < max;
		}
		const Interval i;
		const Int max;
		Int n;
	};

	lib_interval_tree::interval_tree_t<Int> tree_;

};
/****
DIAMOND protein aligner
Copyright (C) 2020 Max Planck Society for the Advancement of Science e.V.

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
#include "../basic/config.h"
#include "../basic/sequence.h"

template<typename _t>
void culling(std::vector<_t> &targets, int source_query_len, const char* query_title, const sequence& query_seq) {
	if (config.min_id > 0 || config.query_cover > 0 || config.subject_cover > 0 || config.no_self_hits)
		for (_t& match : targets)
			match.apply_filters(source_query_len, query_title, query_seq);

	std::sort(targets.begin(), targets.end());
	if (targets.empty() || targets.front().filter_score == 0) {
		targets.clear();
		return;
	}

	typename std::vector<_t>::iterator i = targets.begin();
	if (config.toppercent < 100.0) {
		const int cutoff = std::max(top_cutoff_score(targets.front().filter_score), 1);
		while (i < targets.end() && i->filter_score >= cutoff)
			++i;
	}
	else {
		i += std::min((size_t)config.max_alignments, targets.size());
		while (--i > targets.begin() && i->filter_score == 0);
		++i;
	}
	targets.erase(i, targets.end());
}
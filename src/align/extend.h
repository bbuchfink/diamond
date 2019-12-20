/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen

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

#ifndef EXTEND_H_
#define EXTEND_H_

#include <list>
#include <algorithm>
#include "../basic/parameters.h"
#include "../search/trace_pt_buffer.h"
#include "../data/metadata.h"
#include "../basic/match.h"

namespace Extension {

struct Match {
	Match():
		filter_score(0)
	{}
	void add_hit(std::list<Hsp> &list, std::list<Hsp>::iterator it) {
		hsp.splice(hsp.end(), list, it);
	}
	bool operator<(const Match &m) const {
		return filter_score > m.filter_score;
	}
	void inner_culling(int source_query_len);
	void apply_filters(int source_query_len, const char *query_title);
	int filter_score;
	std::list<Hsp> hsp;
};
	
std::vector<Match> extend(const Parameters &params, size_t query_id, Trace_pt_list::iterator begin, Trace_pt_list::iterator end, const Metadata &metadata);

}

#endif
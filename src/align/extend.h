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
#include <list>
#include <algorithm>
#include <float.h>
#include "../basic/parameters.h"
#include "../search/trace_pt_buffer.h"
#include "../data/metadata.h"
#include "../basic/match.h"
#include "../basic/statistics.h"
#include "../util/text_buffer.h"

namespace Extension {

struct Match {
	Match(size_t target_block_id, int ungapped_score, int filter_score = 0, double filter_evalue = DBL_MAX):
		target_block_id(target_block_id),
		filter_score(filter_score),
		filter_evalue(filter_evalue),
		ungapped_score(ungapped_score)
	{}
	void add_hit(std::list<Hsp> &list, std::list<Hsp>::iterator it) {
		hsp.splice(hsp.end(), list, it);
	}
	static bool cmp_evalue(const Match& m, const Match& n) {
		return m.filter_evalue < n.filter_evalue || (m.filter_evalue == n.filter_evalue && cmp_score(m, n));
	}
	static bool cmp_score(const Match& m, const Match& n) {
		return m.filter_score > n.filter_score || (m.filter_score == n.filter_score && m.target_block_id < n.target_block_id);
	}
	Match(size_t target_block_id, std::array<std::list<Hsp>, MAX_CONTEXT> &hsp, int ungapped_score);
	void inner_culling(int source_query_len);
	void max_hsp_culling();
	void apply_filters(int source_query_len, const char *query_title, const sequence& query_seq);
	size_t target_block_id;
	int filter_score;
	double filter_evalue;
	int ungapped_score;
	std::list<Hsp> hsp;
};

std::vector<Match> extend(const Parameters &params, size_t query_id, hit* begin, hit* end, const Metadata &metadata, Statistics &stat, int flags);
TextBuffer* generate_output(vector<Match> &targets, size_t query_block_id, Statistics &stat, const Metadata &metadata, const Parameters &parameters);
TextBuffer* generate_intermediate_output(const vector<Match> &targets, size_t query_block_id);

/*inline int raw_score_cutoff(size_t query_len) {
	return score_matrix.rawscore(config.min_bit_score == 0 ? score_matrix.bitscore(config.max_evalue, (unsigned)query_len) : config.min_bit_score);
}*/

}

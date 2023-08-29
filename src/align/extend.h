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
#include "../basic/match.h"
#include "../basic/statistics.h"
#include "../util/text_buffer.h"
#include "../run/config.h"
#include "../dp/flags.h"
#include "../stats/cbs.h"
#include "../output/def.h"

namespace Extension {

extern const std::map<Sensitivity, Mode> default_ext_mode;

struct Match {
	Match(const BlockId target_block_id, const Sequence& seq, const ::Stats::TargetMatrix& matrix, Score ungapped_score, Score filter_score = 0, double filter_evalue = DBL_MAX):
		target_block_id(target_block_id),
		seq(seq),
		matrix(matrix),
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
	Match(BlockId target_block_id, const Sequence& seq, const ::Stats::TargetMatrix& matrix, std::array<std::list<Hsp>, MAX_CONTEXT> &hsp, int ungapped_score);
	static Match self_match(BlockId query_id, Sequence query_seq);
	void inner_culling();
	void max_hsp_culling();
	void apply_filters(int source_query_len, const char *query_title, const Sequence& query_seq, const double query_self_aln_score, const Block& targets, const OutputFormat* output_format);
	BlockId target_block_id;
	Sequence seq;
	::Stats::TargetMatrix matrix;
	int filter_score;
	double filter_evalue;
	int ungapped_score;
	std::list<Hsp> hsp;
};

std::pair<std::vector<Match>, Stats> extend(BlockId query_id, Search::Hit* begin, Search::Hit* end, const Search::Config &cfg, Statistics &stat, DP::Flags flags);
TextBuffer* generate_output(std::vector<Match> &targets, const Extension::Stats& stats, BlockId query_block_id, Statistics &stat, const Search::Config& cfg);
TextBuffer* generate_intermediate_output(const std::vector<Match> &targets, BlockId query_block_id, const Search::Config& cfg);

}

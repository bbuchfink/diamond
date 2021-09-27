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

namespace Extension {

enum class Mode {
	BANDED_FAST, BANDED_SLOW, FULL, GLOBAL
};

extern const std::map<Sensitivity, Mode> default_ext_mode;

struct Match {
	Match(size_t target_block_id, const Sequence& seq, const Stats::TargetMatrix& matrix, int ungapped_score, int filter_score = 0, double filter_evalue = DBL_MAX):
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
	Match(size_t target_block_id, const Sequence& seq, const Stats::TargetMatrix& matrix, std::array<std::list<Hsp>, MAX_CONTEXT> &hsp, int ungapped_score);
	void inner_culling();
	void max_hsp_culling();
	void apply_filters(int source_query_len, const char *query_title, const Sequence& query_seq, const Block& targets);
	size_t target_block_id;
	Sequence seq;
	Stats::TargetMatrix matrix;
	int filter_score;
	double filter_evalue;
	int ungapped_score;
	std::list<Hsp> hsp;
};

std::vector<Match> extend(size_t query_id, Search::Hit* begin, Search::Hit* end, const Search::Config &cfg, Statistics &stat, DP::Flags flags);
TextBuffer* generate_output(vector<Match> &targets, size_t query_block_id, Statistics &stat, const Search::Config& cfg);
TextBuffer* generate_intermediate_output(const vector<Match> &targets, size_t query_block_id, const Search::Config& cfg);

/*inline int raw_score_cutoff(size_t query_len) {
	return score_matrix.rawscore(config.min_bit_score == 0 ? score_matrix.bitscore(config.max_evalue, (unsigned)query_len) : config.min_bit_score);
}*/

}

template<>
struct EnumTraits<Extension::Mode> {
	static const SEMap<Extension::Mode> from_string;
};

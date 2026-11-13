/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once
#include <cstdint>
#include <list>
#include <array>
#include <vector>
#include "basic/match.h"
#include "basic/config.h"
#include "basic/statistics.h"
#include "util/text_buffer.h"
#include "run/config.h"
#include "dp/flags.h"
#include "stats/cbs.h"
#include "stats/hauser_correction.h"
#include "basic/const.h"
#include "util/sequence/seqindex.h"

namespace Extension {

extern const std::map<Sensitivity, Mode> default_ext_mode;

struct Query {
	Query(BlockId block_id, Statistics& stats, const Search::Config& cfg, std::pmr::memory_resource& pool);

	const int8_t* composition_bias(int context) const {
		return hauser_correction.empty() ? nullptr : hauser_correction[context].int8.data();
	}

	BlockId block_id;
	std::vector<Sequence> sequence;
	std::vector<HauserCorrection> hauser_correction;
	::Stats::Composition composition{};
	const char* title;
	int source_length;
	int true_aa_length;
	double self_alignment_score;
	std::unique_ptr<Seqindex> seqindex;
	Score ungapped_cutoff;
};

struct Match {
	Match(const BlockId target_block_id, const Sequence& seq, std::unique_ptr<::Stats::TargetMatrix>&& matrix, Score ungapped_score, Score filter_score = 0, double filter_evalue = DBL_MAX):
		target_block_id(target_block_id),
		seq(seq),
		matrix(std::move(matrix)),
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
		if (config.float_scores)
			return m.filter_bitscore() > n.filter_bitscore()
				|| (m.filter_bitscore() == n.filter_bitscore() && m.target_block_id < n.target_block_id);
		return m.filter_score > n.filter_score || (m.filter_score == n.filter_score && m.target_block_id < n.target_block_id);
	}
	double filter_bitscore() const {
		return config.float_scores ? (hsp.empty() ? 0.0 : hsp.front().bit_score) : score_matrix.bitscore(filter_score);
	}
	Match(BlockId target_block_id, const Sequence& seq, std::unique_ptr<::Stats::TargetMatrix>&& matrix, std::array<std::list<Hsp>, MAX_CONTEXT> &hsp, int ungapped_score);
	static Match self_match(BlockId query_id, Sequence query_seq);
	void inner_culling();
	void max_hsp_culling();
	void apply_filters(const Query& query, const Block& targets, const OutputFormat* output_format);
	BlockId target_block_id;
	Sequence seq;
	std::unique_ptr<::Stats::TargetMatrix> matrix;
	int filter_score;
	double filter_evalue;
	int ungapped_score;
	std::list<Hsp> hsp;
};

std::vector<Match> extend(BlockId query_id, Search::Hit* begin, Search::Hit* end, const Search::Config &cfg, Statistics &stat, DP::Flags flags, std::pmr::memory_resource& pool);
std::vector<Match> extend_floating(const Query& query, const Search::Config& cfg);
TextBuffer* generate_output(std::vector<Match> &targets, BlockId query_block_id, Statistics &stat, const Search::Config& cfg);
TextBuffer* generate_intermediate_output(const std::vector<Match> &targets, BlockId query_block_id, const Search::Config& cfg);

}

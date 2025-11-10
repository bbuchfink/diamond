/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <list>
#include <array>
#include "basic/match.h"
#include "basic/statistics.h"
#include "util/text_buffer.h"
#include "run/config.h"
#include "dp/flags.h"
#include "stats/cbs.h"
#include "basic/const.h"

namespace Extension {

extern const std::map<Sensitivity, Mode> default_ext_mode;

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
		return m.filter_score > n.filter_score || (m.filter_score == n.filter_score && m.target_block_id < n.target_block_id);
	}
	Match(BlockId target_block_id, const Sequence& seq, std::unique_ptr<::Stats::TargetMatrix>&& matrix, std::array<std::list<Hsp>, MAX_CONTEXT> &hsp, int ungapped_score);
	static Match self_match(BlockId query_id, Sequence query_seq);
	void inner_culling();
	void max_hsp_culling();
	void apply_filters(int source_query_len, const char *query_title, const Sequence& query_seq, const double query_self_aln_score, const Block& targets, const OutputFormat* output_format);
	BlockId target_block_id;
	Sequence seq;
	std::unique_ptr<::Stats::TargetMatrix> matrix;
	int filter_score;
	double filter_evalue;
	int ungapped_score;
	std::list<Hsp> hsp;
};

std::vector<Match> extend(BlockId query_id, Search::Hit* begin, Search::Hit* end, const Search::Config &cfg, Statistics &stat, DP::Flags flags, std::pmr::monotonic_buffer_resource& pool);
TextBuffer* generate_output(std::vector<Match> &targets, BlockId query_block_id, Statistics &stat, const Search::Config& cfg);
TextBuffer* generate_intermediate_output(const std::vector<Match> &targets, BlockId query_block_id, const Search::Config& cfg);

}
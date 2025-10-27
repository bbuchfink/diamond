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

#include "output/output.h"
#include "basic/statistics.h"
#include "util/data_structures/bit_vector.h"
#include "../extend.h"
#include "../target.h"

namespace Extension { namespace GlobalRanking {

struct QueryList {

	struct Target {
		uint32_t database_id;
		uint16_t score;
	};

	uint32_t query_block_id, last_query_block_id;

	std::vector<Target> targets;
};

std::vector<Extension::Match> ranking_list(size_t query_id, std::vector<TargetScore>::iterator begin, std::vector<TargetScore>::iterator end, std::vector<uint32_t>::const_iterator target_block_ids, const FlatArray<SeedHit>& seed_hits, const Search::Config& cfg);
void write_merged_query_list(const IntermediateRecord& r, TextBuffer& out, BitVector& ranking_db_filter, Statistics& stat);
size_t write_merged_query_list_intro(uint32_t query_id, TextBuffer& buf);
void finish_merged_query_list(TextBuffer& buf, size_t seek_pos);
void extend(SequenceFile& db, TempFile& merged_query_list, BitVector& ranking_db_filter, Search::Config& cfg, Consumer& master_out);
void extend(Search::Config& cfg, Consumer& out);
QueryList fetch_query_targets(InputFile& query_list, uint32_t& next_query);
void update_table(Search::Config& cfg);

struct Hit {
	Hit():
		oid(),
		score(0)
	{}
	Hit(uint32_t oid, uint16_t score, unsigned context):
		oid(oid),
		score(score),
		context((uint8_t)context)
	{}
	Hit(ptrdiff_t target_id):
		oid((uint32_t)target_id),
		score()
	{}
	uint32_t oid;
	uint16_t score;
	uint8_t context;
	bool operator<(const Hit& x) const {
		return score > x.score || (score == x.score && oid < x.oid);
	}
	struct Target {
		uint32_t operator()(const Hit& h) const {
			return h.oid;
		}
	};
	struct CmpOidScore {
		bool operator()(const Hit&x, const Hit& y) const {
			return x.oid < y.oid || (x.oid == y.oid && x.score > y.score);
		}
	};
	struct CmpOid {
		bool operator()(const Hit&x, const Hit& y) const {
			return x.oid == y.oid;
		}
	};
};

}}
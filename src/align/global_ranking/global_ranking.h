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

#include "../../output/output.h"
#include "../../basic/statistics.h"
#include "../../util/data_structures/bit_vector.h"
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
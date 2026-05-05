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
#include "util/data_structures/double_array.h"
#include "data/flags.h"
#include "data/sequence_set.h"
#include "basic/seed.h"
#include "data/block/block.h"

namespace Search {

struct KmerRanking {

	KmerRanking(const SequenceSet& queries, SeedPartition seedp_count, DoubleArray<PackedLocId> *query_seed_hits, DoubleArray<PackedLocId> *ref_seed_hits);
	KmerRanking(const SequenceSet& queries, SeedPartition seedp_count, DoubleArray<PackedLoc>* query_seed_hits, DoubleArray<PackedLoc>* ref_seed_hits);
	KmerRanking(const SequenceSet& queries);

	int highest_ranking(const PackedLocId* begin, const PackedLocId* end, const Block* block) {
		ptrdiff_t r = 0;
		float rank = rank_[begin->block_id];
		uint64_t id;
		if (block)
			id = std::stoll(block->ids()[begin->block_id]);
		
		for (const PackedLocId* i = begin + 1; i < end; ++i) {
			uint64_t this_id;
			if(block) this_id = std::stoll(block->ids()[i->block_id]);
			if (rank_[i->block_id] > rank || (rank_[i->block_id] == rank && block && this_id < id)) {
				rank = rank_[i->block_id];
				id = this_id;
				r = i - begin;
			}
		}
		return (int)r;
	}

private:

	std::vector<float> rank_;

};

}
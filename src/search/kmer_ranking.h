/****
Copyright (C) 2012-2026 Benjamin J. Buchfink <buchfink@gmail.com>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

	http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
****/
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include "util/data_structures/double_array.h"
#include "data/flags.h"
#include "data/sequence_set.h"
#include "basic/seed.h"

namespace Search {

struct KmerRanking {

	KmerRanking(const SequenceSet& queries, SeedPartition seedp_count, DoubleArray<PackedLocId> *query_seed_hits, DoubleArray<PackedLocId> *ref_seed_hits);
	KmerRanking(const SequenceSet& queries, SeedPartition seedp_count, DoubleArray<PackedLoc>* query_seed_hits, DoubleArray<PackedLoc>* ref_seed_hits);
	KmerRanking(const SequenceSet& queries);

	int highest_ranking(const PackedLocId* begin, const PackedLocId* end) {
		ptrdiff_t r = 0;
		float rank = rank_[begin->block_id];
		for (const PackedLocId* i = begin + 1; i < end; ++i)
			if (rank_[i->block_id] > rank) {
				rank = rank_[i->block_id];
				r = i - begin;
			}
		return (int)r;
	}

private:

	std::vector<float> rank_;

};

}
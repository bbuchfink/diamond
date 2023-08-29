#pragma once
#include "../util/data_structures/double_array.h"
#include "../data/flags.h"
#include "../data/sequence_set.h"

namespace Search {

struct KmerRanking {

	KmerRanking(const SequenceSet& queries, DoubleArray<SeedLoc> *query_seed_hits, DoubleArray<SeedLoc> *ref_seed_hits);
	KmerRanking(const SequenceSet& queries);

	int highest_ranking(const SeedLoc* begin, const SeedLoc* end) {
		ptrdiff_t r = 0;
#ifdef KEEP_TARGET_ID
		float rank = rank_[begin->block_id];
		for (const SeedLoc* i = begin + 1; i < end; ++i)
			if (rank_[i->block_id] > rank) {
				rank = rank_[i->block_id];
				r = i - begin;
			}
#endif
		return (int)r;
	}

private:

	std::vector<float> rank_;

};

}
#pragma once
#include "../util/data_structures/double_array.h"
#include "../data/flags.h"
#include "../data/sequence_set.h"

namespace Search {

struct KmerRanking {

	KmerRanking(const SequenceSet& queries, DoubleArray<PackedLocId> *query_seed_hits, DoubleArray<PackedLocId> *ref_seed_hits);
	KmerRanking(const SequenceSet& queries, DoubleArray<PackedLoc>* query_seed_hits, DoubleArray<PackedLoc>* ref_seed_hits);
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
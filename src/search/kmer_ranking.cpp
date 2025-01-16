#include <atomic>
#include <thread>
#include "kmer_ranking.h"
#include "util/algo/join_result.h"
#include "basic/config.h"

using std::pair;
using std::vector;
using std::atomic_int;
using std::thread;
using std::greater;
using std::atomic;
using std::runtime_error;

namespace Search {

KmerRanking::KmerRanking(const SequenceSet& queries, SeedPartition seedp_count, DoubleArray<PackedLoc>* query_seed_hits, DoubleArray<PackedLoc>* ref_seed_hits) {
	throw runtime_error("Unsupported");
}

KmerRanking::KmerRanking(const SequenceSet& queries, SeedPartition seedp_count, DoubleArray<PackedLocId> *query_seed_hits, DoubleArray<PackedLocId> *ref_seed_hits) {
	atomic<SeedPartition> rel_seedp(0);
	const auto query_count = queries.size();
	atomic<float>* counts = new atomic<float>[query_count];
	
	auto worker = [&] {
		SeedPartition p;
		while ((p = rel_seedp.fetch_add(1, std::memory_order_relaxed)) < seedp_count) {
			for (auto it = JoinIterator<PackedLocId>(query_seed_hits[p].begin(), ref_seed_hits[p].begin()); it;) {
				const Range<PackedLocId*> query_hits = *it.r;
				for (PackedLocId s : query_hits) {
					float value = counts[s.block_id].load();
					while (!counts[s.block_id].compare_exchange_weak(value, value + (float)sqrt(it.s->size())));
					//counts[s.block_id] += it.s->size();
				}
				++it;
			}
		}
		
	};

	vector<thread> threads;
	for (int i = 0; i < config.threads_; ++i)
		threads.emplace_back(worker);
	for (auto& i : threads)
		i.join();
	
	rank_.reserve(query_count);
	for (BlockId i = 0; i < query_count; ++i)
		rank_.push_back(counts[i]);
	delete[] counts;
}

KmerRanking::KmerRanking(const SequenceSet& queries) {
	rank_.reserve(queries.size());
	for (BlockId i = 0; i < queries.size(); ++i)
		rank_.push_back((float)queries[i].length());
}

}
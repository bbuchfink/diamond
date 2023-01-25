#include <atomic>
#include <thread>
#include "kmer_ranking.h"
#include "../util/algo/join_result.h"
#include "../basic/config.h"

using std::pair;
using std::vector;
using std::atomic_int;
using std::thread;
using std::greater;
using std::atomic;
using std::runtime_error;

namespace Search {

KmerRanking::KmerRanking(const SequenceSet& queries, DoubleArray<SeedLoc> *query_seed_hits, DoubleArray<SeedLoc> *ref_seed_hits) {
#ifdef KEEP_TARGET_ID
	atomic_int seedp(0);
	const auto query_count = queries.size();
	atomic<float>* counts = new atomic<float>[query_count];
	
	auto worker = [&] {
		int p;
		while ((p = seedp++) < Const::seedp) {
			for (auto it = JoinIterator<SeedLoc>(query_seed_hits[p].begin(), ref_seed_hits[p].begin()); it;) {
				const Range<SeedLoc*> query_hits = *it.r;
				for (SeedLoc s : query_hits) {
					float value = counts[s.block_id].load();
					while (!counts[s.block_id].compare_exchange_weak(value, value + (float)sqrt(it.s->size())));
					//counts[s.block_id] += it.s->size();
				}
				++it;
			}
		}
		
	};

	vector<thread> threads;
	for (size_t i = 0; i < config.threads_; ++i)
		threads.emplace_back(worker);
	for (auto& i : threads)
		i.join();
	
	rank_.reserve(query_count);
	for (BlockId i = 0; i < query_count; ++i)
		rank_.push_back(counts[i]);
	delete[] counts;
#else
	throw runtime_error("Option only available when compiled with KEEP_TARGET_ID=ON");
#endif
}

KmerRanking::KmerRanking(const SequenceSet& queries) {
	rank_.reserve(queries.size());
	for (BlockId i = 0; i < queries.size(); ++i)
		rank_.push_back((float)queries[i].length());
}

}
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

#include <numeric>
#include <atomic>
#include "frequent_seeds.h"
#include "util/parallel/thread_pool.h"
#include "util/util.h"
#include "basic/config.h"
#include "basic/seed.h"
#include "sequence_set.h"
#include "block/block.h"
#include "util/algo/join_result.h"

using std::endl;
using std::atomic;
using std::vector;

const double FrequentSeeds::hash_table_factor = 1.3;
FrequentSeeds frequent_seeds;

template<typename SeedLoc>
static void compute_sd(atomic<SeedPartition> *seedp, SeedPartition partition_count, DoubleArray<SeedLoc> *query_seed_hits, DoubleArray<SeedLoc> *ref_seed_hits, vector<Sd> *ref_out, vector<Sd> *query_out)
{
	SeedPartition p;
	while ((p = seedp->fetch_add(1, std::memory_order_relaxed)) < partition_count) {
		Sd ref_sd, query_sd;
		for (auto it = JoinIterator<SeedLoc>(query_seed_hits[p].begin(), ref_seed_hits[p].begin()); it; ++it) {
			query_sd.add((double)it.r->size());
			ref_sd.add((double)it.s->size());
		}
		(*ref_out)[p] = ref_sd;
		(*query_out)[p] = query_sd;
	}
}

template<typename SeedLoc>
void FrequentSeeds::build_worker(
	size_t rel_seedp,
	size_t thread_id,
	DoubleArray<SeedLoc> *query_seed_hits,
	DoubleArray<SeedLoc> *ref_seed_hits,
	const SeedPartitionRange *range,
	unsigned sid,
	unsigned ref_max_n,
	unsigned query_max_n,
	vector<unsigned> *counts,
	Search::Config* cfg)
{
	SequenceSet& query_seqs = cfg->query->seqs();
	
	vector<uint32_t> buf;
	size_t n = 0;
	for (auto it = JoinIterator<SeedLoc>(query_seed_hits[rel_seedp].begin(), ref_seed_hits[rel_seedp].begin()); it;) {
		if (it.s->size() > ref_max_n || it.r->size() > query_max_n) {
			n += (unsigned)it.s->size();
			//Packed_seed s;
			//shapes[sid].set_seed(s, query_seqs::get().data(*it.r->begin()));
			//buf.push_back(seed_partition_offset(s));

			Range<SeedLoc*> query_hits = *it.r;
			for (SeedLoc* i = query_hits.begin(); i < query_hits.end(); ++i) {
				Letter* p = query_seqs.data(*i);
				*p |= SEED_MASK;
			}
			
			it.erase();

		}
		else
			++it;
	}

	(*counts)[rel_seedp] = (unsigned)n;
}

template<typename SeedLoc>
void FrequentSeeds::build(unsigned sid, const SeedPartitionRange &range, DoubleArray<SeedLoc> *query_seed_hits, DoubleArray<SeedLoc> *ref_seed_hits, Search::Config& cfg)
{
	vector<Sd> ref_sds(range.size()), query_sds(range.size());
	atomic<unsigned> rel_seedp(0);
	vector<std::thread> threads;
	for (int i = 0; i < config.threads_; ++i)
		threads.emplace_back(compute_sd<SeedLoc>, &rel_seedp, range.size(), query_seed_hits, ref_seed_hits, &ref_sds, &query_sds);
	for (auto &t : threads)
		t.join();

	Sd ref_sd(ref_sds), query_sd(query_sds);
	const unsigned ref_max_n = (unsigned)(ref_sd.mean() + cfg.freq_sd*ref_sd.sd()), query_max_n = (unsigned)(query_sd.mean() + cfg.freq_sd*query_sd.sd());
	log_stream << "Seed frequency mean (reference) = " << ref_sd.mean() << ", SD = " << ref_sd.sd() << endl;
	log_stream << "Seed frequency mean (query) = " << query_sd.mean() << ", SD = " << query_sd.sd() << endl;
	log_stream << "Seed frequency cap query: " << query_max_n << ", reference: " << ref_max_n << endl;
	vector<unsigned> counts(range.size());
	Util::Parallel::scheduled_thread_pool_auto(config.threads_, range.size(), build_worker<SeedLoc>, query_seed_hits, ref_seed_hits, &range, sid, ref_max_n, query_max_n, &counts, &cfg);
	log_stream << "Masked positions = " << std::accumulate(counts.begin(), counts.end(), 0) << std::endl;
}

template void FrequentSeeds::build(unsigned, const SeedPartitionRange&, DoubleArray<PackedLoc>*, DoubleArray<PackedLoc>*, Search::Config&);
template void FrequentSeeds::build(unsigned, const SeedPartitionRange&, DoubleArray<PackedLocId>*, DoubleArray<PackedLocId>*, Search::Config&);

void FrequentSeeds::clear_masking(SequenceSet& seqs) {
	for (BlockId i = 0; i < seqs.size(); ++i) {
		const size_t len = seqs.length(i);
		Letter* p = seqs.ptr(i), *end = p + len;
		for (; p < end; ++p)
			*p = letter_mask(*p);
	}
}
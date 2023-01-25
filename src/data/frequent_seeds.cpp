/****
DIAMOND protein aligner
Copyright (C) 2016-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
						
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

#include <numeric>
#include <utility>
#include <atomic>
#include "frequent_seeds.h"
#include "queries.h"
#include "../util/parallel/thread_pool.h"
#include "../util/util.h"
#include "../basic/config.h"

using std::endl;
using std::atomic;
using std::vector;

const double Frequent_seeds::hash_table_factor = 1.3;
Frequent_seeds frequent_seeds;

static void compute_sd(atomic<unsigned> *seedp, DoubleArray<SeedLoc> *query_seed_hits, DoubleArray<SeedLoc> *ref_seed_hits, vector<Sd> *ref_out, vector<Sd> *query_out)
{
	int p;
	while ((p = (*seedp)++) < current_range.end()) {
		Sd ref_sd, query_sd;
		for (auto it = JoinIterator<SeedLoc>(query_seed_hits[p].begin(), ref_seed_hits[p].begin()); it; ++it) {
			query_sd.add((double)it.r->size());
			ref_sd.add((double)it.s->size());
		}
		(*ref_out)[p - current_range.begin()] = ref_sd;
		(*query_out)[p - current_range.begin()] = query_sd;
	}
}

void Frequent_seeds::build_worker(
	size_t seedp,
	size_t thread_id,
	DoubleArray<SeedLoc> *query_seed_hits,
	DoubleArray<SeedLoc> *ref_seed_hits,
	const SeedPartitionRange *range,
	unsigned sid,
	unsigned ref_max_n,
	unsigned query_max_n,
	vector<unsigned> *counts,
	Search::Config* cfg) {
	SequenceSet& query_seqs = cfg->query->seqs();
	if (!range->contains((unsigned)seedp))
		return;

	vector<uint32_t> buf;
	size_t n = 0;
	for (auto it = JoinIterator<SeedLoc>(query_seed_hits[seedp].begin(), ref_seed_hits[seedp].begin()); it;) {
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

	(*counts)[seedp] = (unsigned)n;
}

void Frequent_seeds::build(unsigned sid, const SeedPartitionRange &range, DoubleArray<SeedLoc> *query_seed_hits, DoubleArray<SeedLoc> *ref_seed_hits, Search::Config& cfg)
{
	vector<Sd> ref_sds(range.size()), query_sds(range.size());
	atomic<unsigned> seedp(range.begin());
	vector<std::thread> threads;
	for (int i = 0; i < config.threads_; ++i)
		threads.emplace_back(compute_sd, &seedp, query_seed_hits, ref_seed_hits, &ref_sds, &query_sds);
	for (auto &t : threads)
		t.join();

	Sd ref_sd(ref_sds), query_sd(query_sds);
	const unsigned ref_max_n = (unsigned)(ref_sd.mean() + cfg.freq_sd*ref_sd.sd()), query_max_n = (unsigned)(query_sd.mean() + cfg.freq_sd*query_sd.sd());
	log_stream << "Seed frequency mean (reference) = " << ref_sd.mean() << ", SD = " << ref_sd.sd() << endl;
	log_stream << "Seed frequency mean (query) = " << query_sd.mean() << ", SD = " << query_sd.sd() << endl;
	log_stream << "Seed frequency cap query: " << query_max_n << ", reference: " << ref_max_n << endl;
	vector<unsigned> counts(Const::seedp);
	Util::Parallel::scheduled_thread_pool_auto(config.threads_, Const::seedp, build_worker, query_seed_hits, ref_seed_hits, &range, sid, ref_max_n, query_max_n, &counts, &cfg);
	log_stream << "Masked positions = " << std::accumulate(counts.begin(), counts.end(), 0) << std::endl;
}

void Frequent_seeds::clear_masking(SequenceSet& seqs) {
	for (BlockId i = 0; i < seqs.size(); ++i) {
		const size_t len = seqs.length(i);
		Letter* p = seqs.ptr(i), *end = p + len;
		for (; p < end; ++p)
			*p = letter_mask(*p);
	}
}
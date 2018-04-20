/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#include "search.h"
#include "../util/algo/hash_join.h"
#include "../util/algo/radix_sort.h"
#include "../data/reference.h"
#include "../data/seed_array.h"
#include "../data/queries.h"
#include "../data/frequent_seeds.h"
#include "trace_pt_buffer.h"
#include "align_range.h"

void seed_join_worker(SeedArray *query_seeds, SeedArray *ref_seeds, Atomic<unsigned> *seedp, const SeedPartitionRange *seedp_range, typename vector<JoinResult<SeedArray::Entry> >::iterator seed_hits)
{
	unsigned p;
	MemoryPool tmp_pool(false);
	while ((p = (*seedp)++) < seedp_range->end()) {
		hash_join(
			Relation<SeedArray::Entry>(query_seeds->begin(p), query_seeds->size(p)),
			Relation<SeedArray::Entry>(ref_seeds->begin(p), ref_seeds->size(p)),
			*(seed_hits + (p - seedp_range->begin())),
			tmp_pool,
			24);
		MemoryPool::global().free(query_seeds->begin(p));
		MemoryPool::global().free(ref_seeds->begin(p));
	}
}

void search_worker(Atomic<unsigned> *seedp, const SeedPartitionRange *seedp_range, unsigned shape, size_t thread_id, vector<JoinResult<SeedArray::Entry> >::iterator seed_hits)
{
	Trace_pt_buffer::Iterator* out = new Trace_pt_buffer::Iterator(*Trace_pt_buffer::instance, thread_id);
	Statistics stats;
	Seed_filter seed_filter(stats, *out, shape);
	unsigned p;
	while ((p = (*seedp)++) < seedp_range->end())
		for (JoinResult<SeedArray::Entry>::Iterator it = seed_hits[p - seedp_range->begin()].begin(); it.good(); ++it)
			if (it.s[0] != 0)
				seed_filter.run(it.r.data(), it.r.count(), it.s.data(), it.s.count());
	delete out;
	statistics += stats;
}

void search_shape(unsigned sid, unsigned query_block)
{
	::partition<unsigned> p(Const::seedp, config.lowmem);
	for (unsigned chunk = 0; chunk < p.parts; ++chunk) {
		message_stream << "Processing query block " << query_block << ", reference block " << current_ref_block << ", shape " << sid << ", index chunk " << chunk << '.' << endl;
		const SeedPartitionRange range(p.getMin(chunk), p.getMax(chunk));
		current_range = range;

		task_timer timer("Building reference seed array", true);
		SeedArray *ref_idx;
		if (config.algo == Config::query_indexed)
			ref_idx = new SeedArray(*ref_seqs::data_, sid, ref_hst.get(sid), range, ref_hst.partition(), query_seeds);
		else if (query_seeds_hashed != 0)
			ref_idx = new SeedArray(*ref_seqs::data_, sid, ref_hst.get(sid), range, ref_hst.partition(), query_seeds_hashed);
		else
			ref_idx = new SeedArray(*ref_seqs::data_, sid, ref_hst.get(sid), range, ref_hst.partition(), &no_filter);

		timer.go("Building query seed array");
		SeedArray query_idx(*query_seqs::data_, sid, query_hst.get(sid), range, query_hst.partition(), &no_filter);

		timer.go("Computing hash join");
		//MemoryPool::init((query_seeds.limits_.back() + ref_seeds.limits_.back()) * 5);
		Atomic<unsigned> seedp(range.begin());
		Thread_pool threads;
		vector<JoinResult<SeedArray::Entry> > seed_hits(range.size());
		for (size_t i = 0; i < config.threads_; ++i)
			threads.push_back(launch_thread(seed_join_worker, &query_idx, ref_idx, &seedp, &range, seed_hits.begin()));
		threads.join_all();
		delete ref_idx;

		timer.go("Building seed filter");
		frequent_seeds.build(sid, range, seed_hits.begin());

		timer.go("Searching alignments");
		seedp = range.begin();
		for (size_t i = 0; i < config.threads_; ++i)
			threads.push_back(launch_thread(search_worker, &seedp, &range, sid, i, seed_hits.begin()));
		threads.join_all();
	}
}
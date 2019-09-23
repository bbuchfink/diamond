/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

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

#include <thread>
#include <utility>
#include "search.h"
#include "../util/algo/hash_join.h"
#include "../util/algo/radix_sort.h"
#include "../data/reference.h"
#include "../data/seed_array.h"
#include "../data/queries.h"
#include "../data/frequent_seeds.h"
#include "trace_pt_buffer.h"
#include "align_range.h"
#include "../util/data_structures/double_array.h"

using namespace std;

void seed_join_worker(
	SeedArray *query_seeds,
	SeedArray *ref_seeds,
	Atomic<unsigned> *seedp,
	const SeedPartitionRange *seedp_range,
	DoubleArray<SeedArray::_pos> *query_seed_hits,
	DoubleArray<SeedArray::_pos> *ref_seeds_hits)
{
	unsigned p;
	while ((p = (*seedp)++) < seedp_range->end()) {
		std::pair<DoubleArray<SeedArray::_pos>, DoubleArray<SeedArray::_pos>> join = hash_join(
			Relation<SeedArray::Entry>(query_seeds->begin(p), query_seeds->size(p)),
			Relation<SeedArray::Entry>(ref_seeds->begin(p), ref_seeds->size(p)),
			24);
		query_seed_hits[p] = join.first;
		ref_seeds_hits[p] = join.second;
	}
}

void search_worker(Atomic<unsigned> *seedp, const SeedPartitionRange *seedp_range, unsigned shape, size_t thread_id, DoubleArray<SeedArray::_pos> *query_seed_hits, DoubleArray<SeedArray::_pos> *ref_seed_hits)
{
	Trace_pt_buffer::Iterator* out = new Trace_pt_buffer::Iterator(*Trace_pt_buffer::instance, thread_id);
	Statistics stats;
	Seed_filter seed_filter(stats, *out, shape);
	unsigned p;
	while ((p = (*seedp)++) < seedp_range->end())
		for (auto it = JoinIterator<SeedArray::_pos>(query_seed_hits[p].begin(), ref_seed_hits[p].begin()); it; ++it)
			seed_filter.run(it.r->begin(), it.r->size(), it.s->begin(), it.s->size());
	delete out;
	statistics += stats;
}

void search_shape(unsigned sid, unsigned query_block)
{
	::partition<unsigned> p(Const::seedp, config.lowmem);
	DoubleArray<SeedArray::_pos> query_seed_hits[Const::seedp], ref_seed_hits[Const::seedp];

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
		SeedArray *query_idx = new SeedArray(*query_seqs::data_, sid, query_hst.get(sid), range, query_hst.partition(), &no_filter);

		timer.go("Computing hash join");
		Atomic<unsigned> seedp(range.begin());
		vector<thread> threads;
		for (size_t i = 0; i < config.threads_; ++i)
			threads.emplace_back(seed_join_worker, query_idx, ref_idx, &seedp, &range, query_seed_hits, ref_seed_hits);
		for (auto &t : threads)
			t.join();

		timer.go("Building seed filter");
		frequent_seeds.build(sid, range, query_seed_hits, ref_seed_hits);

		timer.go("Searching alignments");
		seedp = range.begin();
		threads.clear();
		for (size_t i = 0; i < config.threads_; ++i)
			threads.emplace_back(search_worker, &seedp, &range, sid, i, query_seed_hits, ref_seed_hits);
		for (auto &t : threads)
			t.join();

		timer.go("Deallocating buffers");
		delete ref_idx;
		delete query_idx;
	}
}
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

#include <thread>
#include <utility>
#include <atomic>
#include "search.h"
#include "../util/algo/hash_join.h"
#include "../util/algo/radix_sort.h"
#include "../data/reference.h"
#include "../data/seed_array.h"
#include "../data/queries.h"
#include "../data/frequent_seeds.h"
#include "trace_pt_buffer.h"
#include "../util/data_structures/double_array.h"
#include "../util/system/system.h"

using std::vector;
using std::atomic;
using std::endl;
using std::unique_ptr;

Trace_pt_buffer* Trace_pt_buffer::instance;

void seed_join_worker(
	SeedArray *query_seeds,
	SeedArray *ref_seeds,
	atomic<unsigned> *seedp,
	const SeedPartitionRange *seedp_range,
	DoubleArray<SeedArray::_pos> *query_seed_hits,
	DoubleArray<SeedArray::_pos> *ref_seeds_hits)
{
	unsigned p;
	const unsigned bits = query_seeds->key_bits;
	if (bits != ref_seeds->key_bits)
		throw std::runtime_error("Joining seed arrays with different key lengths.");
	while ((p = (*seedp)++) < seedp_range->end()) {
		std::pair<DoubleArray<SeedArray::_pos>, DoubleArray<SeedArray::_pos>> join = hash_join(
			Relation<SeedArray::Entry>(query_seeds->begin(p), query_seeds->size(p)),
			Relation<SeedArray::Entry>(ref_seeds->begin(p), ref_seeds->size(p)),
			bits);
		query_seed_hits[p] = join.first;
		ref_seeds_hits[p] = join.second;
	}
}

void search_worker(atomic<unsigned> *seedp, const SeedPartitionRange *seedp_range, unsigned shape, size_t thread_id, DoubleArray<SeedArray::_pos> *query_seed_hits, DoubleArray<SeedArray::_pos> *ref_seed_hits, const Search::Context *context)
{
#ifdef __APPLE__
	unique_ptr<Search::WorkSet> work_set(new Search::WorkSet{ *context, shape, {},  {*Trace_pt_buffer::instance, thread_id}, {} });
#else
	unique_ptr<Search::WorkSet> work_set(new Search::WorkSet{ *context, shape, {},  {*Trace_pt_buffer::instance, thread_id}, {}, {}, {} });
#endif
	unsigned p;
	while ((p = (*seedp)++) < seedp_range->end())
		for (auto it = JoinIterator<SeedArray::_pos>(query_seed_hits[p].begin(), ref_seed_hits[p].begin()); it; ++it)
			Search::stage1(it.r->begin(), it.r->size(), it.s->begin(), it.s->size(), *work_set);
	statistics += work_set->stats;
}

void search_shape(unsigned sid, unsigned query_block, char *query_buffer, char *ref_buffer, const Parameters &params, const HashedSeedSet* target_seeds)
{
	::partition<unsigned> p(Const::seedp, config.lowmem);
	DoubleArray<SeedArray::_pos> query_seed_hits[Const::seedp], ref_seed_hits[Const::seedp];
	log_rss();

	for (unsigned chunk = 0; chunk < p.parts; ++chunk) {
		message_stream << "Processing query block " << query_block + 1
			<< ", reference block " << (current_ref_block + 1) << "/" << params.ref_blocks
			<< ", shape " << (sid + 1) << "/" << shapes.count();
		if (config.lowmem > 1)
			message_stream << ", index chunk " << chunk + 1 << "/" << config.lowmem;
		message_stream << '.' << endl;
		const SeedPartitionRange range(p.getMin(chunk), p.getMax(chunk));
		current_range = range;

		task_timer timer("Building reference seed array", true);
		SeedArray *ref_idx;
		if (query_seeds_hashed.get())
			ref_idx = new SeedArray(*ref_seqs::data_, sid, ref_hst.get(sid), range, ref_hst.partition(), ref_buffer, query_seeds_hashed.get(), true);
		else
			ref_idx = new SeedArray(*ref_seqs::data_, sid, ref_hst.get(sid), range, ref_hst.partition(), ref_buffer, &no_filter, target_seeds);

		timer.go("Building query seed array");
		SeedArray* query_idx;
		if (target_seeds)
			query_idx = new SeedArray(*query_seqs::data_, sid, range, target_seeds, true);
		else
			query_idx = new SeedArray(*query_seqs::data_, sid, query_hst.get(sid), range, query_hst.partition(), query_buffer, &no_filter, query_seeds_hashed.get());
		timer.finish();

		log_stream << "Indexed query seeds = " << query_idx->size() << '/' << query_seqs::get().letters() << ", reference seeds = " << ref_idx->size() << '/' << ref_seqs::get().letters() << endl;

		timer.go("Computing hash join");
		atomic<unsigned> seedp(range.begin());
		vector<std::thread> threads;
		for (size_t i = 0; i < config.threads_; ++i)
			threads.emplace_back(seed_join_worker, query_idx, ref_idx, &seedp, &range, query_seed_hits, ref_seed_hits);
		for (auto &t : threads)
			t.join();

		timer.go("Building seed filter");
		frequent_seeds.build(sid, range, query_seed_hits, ref_seed_hits);

		Search::Context* context = nullptr;
		const vector<uint32_t> patterns = shapes.patterns(0, sid + 1);
		context = new Search::Context{ {patterns.data(), patterns.data() + patterns.size() - 1 },
			{patterns.data(), patterns.data() + patterns.size() },
			config.ungapped_evalue,
			config.ungapped_evalue_short,
			score_matrix.rawscore(config.short_query_ungapped_bitscore)
		};

		timer.go("Searching alignments");
		seedp = range.begin();
		threads.clear();
		for (size_t i = 0; i < config.threads_; ++i)
			threads.emplace_back(search_worker, &seedp, &range, sid, i, query_seed_hits, ref_seed_hits, context);
		for (auto &t : threads)
			t.join();

		delete ref_idx;
		delete query_idx;
		delete context;
	}
}
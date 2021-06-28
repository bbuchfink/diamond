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
#include "../util/data_structures/double_array.h"
#include "../util/system/system.h"
#include "../util/data_structures/deque.h"
#include "../util/util.h"

using std::vector;
using std::atomic;
using std::endl;
using std::unique_ptr;
using Search::Hit;

void seed_join_worker(
	SeedArray *query_seeds,
	SeedArray *ref_seeds,
	atomic<unsigned> *seedp,
	const SeedPartitionRange *seedp_range,
	DoubleArray<SeedArray::Entry::Value> *query_seed_hits,
	DoubleArray<SeedArray::Entry::Value> *ref_seeds_hits)
{
	unsigned p;
	const unsigned bits = query_seeds->key_bits;
	if (bits != ref_seeds->key_bits)
		throw std::runtime_error("Joining seed arrays with different key lengths.");
	while ((p = (*seedp)++) < seedp_range->end()) {
		std::pair<DoubleArray<SeedArray::Entry::Value>, DoubleArray<SeedArray::Entry::Value>> join = hash_join(
			Relation<SeedArray::Entry>(query_seeds->begin(p), query_seeds->size(p)),
			Relation<SeedArray::Entry>(ref_seeds->begin(p), ref_seeds->size(p)),
			bits);
		query_seed_hits[p] = join.first;
		ref_seeds_hits[p] = join.second;
	}
}

void search_worker(atomic<unsigned> *seedp, const SeedPartitionRange *seedp_range, unsigned shape, size_t thread_id, DoubleArray<SeedArray::Entry::Value> *query_seed_hits, DoubleArray<SeedArray::Entry::Value> *ref_seed_hits, const Search::Context *context, const Search::Config* cfg)
{
	unique_ptr<Writer<Hit>> writer;
	if (config.global_ranking_targets)
		writer.reset(new AsyncWriter<Hit, Search::Config::RankingBuffer::EXPONENT>(*cfg->global_ranking_buffer));
	else
		writer.reset(new AsyncBuffer<Hit>::Iterator(*cfg->seed_hit_buf, thread_id));
#ifdef __APPLE__
	unique_ptr<Search::WorkSet> work_set(new Search::WorkSet{ *context, *cfg, shape, {}, writer.get(), {} });
#else
	unique_ptr<Search::WorkSet> work_set(new Search::WorkSet{ *context, *cfg, shape, {}, writer.get(), {}, {}, {} });
#endif
	unsigned p;
	while ((p = (*seedp)++) < seedp_range->end())
		for (auto it = JoinIterator<SeedArray::Entry::Value>(query_seed_hits[p].begin(), ref_seed_hits[p].begin()); it; ++it)
			Search::stage1(it.r->begin(), it.r->size(), it.s->begin(), it.s->size(), *work_set);
	statistics += work_set->stats;
}

void search_shape(unsigned sid, unsigned query_block, unsigned query_iteration, char *query_buffer, char *ref_buffer, Search::Config& cfg, const HashedSeedSet* target_seeds)
{
	Partition<unsigned> p(Const::seedp, cfg.index_chunks);
	DoubleArray<SeedArray::Entry::Value> query_seed_hits[Const::seedp], ref_seed_hits[Const::seedp];
	log_rss();
	SequenceSet& ref_seqs = cfg.target->seqs(), query_seqs = cfg.query->seqs();
	const Partitioned_histogram& ref_hst = cfg.target->hst(), query_hst = cfg.query->hst();

	for (unsigned chunk = 0; chunk < p.parts; ++chunk) {
		message_stream << "Processing query block " << query_block + 1;
		if (cfg.iterated())
			message_stream << ", query iteration " << query_iteration + 1;
		message_stream << ", reference block " << (current_ref_block + 1) << "/" << cfg.ref_blocks
			<< ", shape " << (sid + 1) << "/" << shapes.count();
		if (cfg.index_chunks > 1)
			message_stream << ", index chunk " << chunk + 1 << "/" << cfg.index_chunks;
		message_stream << '.' << endl;
		const SeedPartitionRange range(p.begin(chunk), p.end(chunk));
		current_range = range;

		task_timer timer("Building reference seed array", true);
		SeedArray *ref_idx;
		if (query_seeds_bitset.get())
			ref_idx = new SeedArray(ref_seqs, sid, ref_hst.get(sid), range, ref_hst.partition(), ref_buffer, query_seeds_bitset.get(), cfg.seed_encoding, nullptr);
		else if (query_seeds_hashed.get())
			ref_idx = new SeedArray(ref_seqs, sid, ref_hst.get(sid), range, ref_hst.partition(), ref_buffer, query_seeds_hashed.get(), cfg.seed_encoding, nullptr);
			//ref_idx = new SeedArray(ref_seqs, sid, range, query_seeds_hashed.get(), true);
		else
			ref_idx = new SeedArray(ref_seqs, sid, ref_hst.get(sid), range, ref_hst.partition(), ref_buffer, &no_filter, cfg.seed_encoding, nullptr);

		timer.go("Building query seed array");
		SeedArray* query_idx;
		if (target_seeds)
			query_idx = new SeedArray(query_seqs, sid, range, target_seeds, cfg.seed_encoding, nullptr);
		else
			query_idx = new SeedArray(query_seqs, sid, query_hst.get(sid), range, query_hst.partition(), query_buffer, &no_filter, cfg.seed_encoding, cfg.query_skip.get());
		timer.finish();

		log_stream << "Indexed query seeds = " << query_idx->size() << '/' << query_seqs.letters() << ", reference seeds = " << ref_idx->size() << '/' << ref_seqs.letters() << endl;

		timer.go("Computing hash join");
		atomic<unsigned> seedp(range.begin());
		vector<std::thread> threads;
		for (size_t i = 0; i < config.threads_; ++i)
			threads.emplace_back(seed_join_worker, query_idx, ref_idx, &seedp, &range, query_seed_hits, ref_seed_hits);
		for (auto &t : threads)
			t.join();

		timer.go("Building seed filter");
		frequent_seeds.build(sid, range, query_seed_hits, ref_seed_hits, cfg);

		Search::Context* context = nullptr;
		const vector<uint32_t> patterns = shapes.patterns(0, sid + 1);
		context = new Search::Context{ {patterns.data(), patterns.data() + patterns.size() - 1 },
			{patterns.data(), patterns.data() + patterns.size() },
			cfg.ungapped_evalue,
			cfg.ungapped_evalue_short,
			score_matrix.rawscore(config.short_query_ungapped_bitscore)
		};

		timer.go("Searching alignments");
		seedp = range.begin();
		threads.clear();
		for (size_t i = 0; i < config.threads_; ++i)
			threads.emplace_back(search_worker, &seedp, &range, sid, i, query_seed_hits, ref_seed_hits, context, &cfg);
		for (auto &t : threads)
			t.join();

		delete ref_idx;
		delete query_idx;
		delete context;
	}
}
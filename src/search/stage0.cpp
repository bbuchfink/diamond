/****
DIAMOND protein aligner
Copyright (C) 2016-2021 Max Planck Society for the Advancement of Science e.V.
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
#include "data/seed_array.h"
#include "data/frequent_seeds.h"
#include "util/data_structures/double_array.h"
#include "util/system/system.h"
#include "util/data_structures/deque.h"
#include "util/async_buffer.h"
#include "basic/seed.h"
#include "seed_complexity.h"
#include "util/algo/partition.h"
#include "data/block/block.h"
#include "basic/shape_config.h"
#include "data/queries.h"
#include "util/algo/radix_cluster.h"
#include "util/algo/hash_join.h"

using std::vector;
using std::atomic;
using std::endl;
using std::unique_ptr;
using Search::Hit;

namespace Search {

template<typename SeedLoc>
static void seed_join_worker(
	SeedArray<SeedLoc> *query_seeds,
	SeedArray<SeedLoc> *ref_seeds,
	atomic<unsigned> *seedp,
	SeedPartition partition_count,
	DoubleArray<SeedLoc> *query_seed_hits,
	DoubleArray<SeedLoc> *ref_seeds_hits)
{
	SeedPartition p;
	const int bits = query_seeds->key_bits;
	if (bits != ref_seeds->key_bits)
		throw std::runtime_error("Joining seed arrays with different key lengths.");
	while ((p = seedp->fetch_add(1, std::memory_order_relaxed)) < partition_count) {
		std::pair<DoubleArray<SeedLoc>, DoubleArray<SeedLoc>> join = hash_join(
			Relation<typename SeedArray<SeedLoc>::Entry>(query_seeds->begin(p), query_seeds->size(p)),
			Relation<typename SeedArray<SeedLoc>::Entry>(ref_seeds->begin(p), ref_seeds->size(p)),
			bits);
		query_seed_hits[p] = join.first;
		ref_seeds_hits[p] = join.second;
	}
}

template<typename SeedLoc>
static void search_worker(atomic<SeedPartition> *seedp, SeedPartition partition_count, unsigned shape, size_t thread_id, DoubleArray<SeedLoc> *query_seed_hits, DoubleArray<SeedLoc> *ref_seed_hits, const Search::Context *context, const Search::Config* cfg)
{
	unique_ptr<Writer<Hit>> writer;
	if (config.global_ranking_targets)
		writer.reset(new AsyncWriter<Hit, Search::Config::RankingBuffer::EXPONENT>(*cfg->global_ranking_buffer));
	else
		writer.reset(new AsyncBuffer<Hit>::Iterator(*cfg->seed_hit_buf, thread_id));
	unique_ptr<Search::WorkSet> work_set(new Search::WorkSet(*context, *cfg, shape, writer.get(), context->kmer_ranking));
	SeedPartition p;
	while ((p = seedp->fetch_add(1, std::memory_order_relaxed)) < partition_count) {
		auto it = JoinIterator<SeedLoc>(query_seed_hits[p].begin(), ref_seed_hits[p].begin());
		run_stage1(it, work_set.get(), cfg);
	}
	statistics += work_set->stats;
}

template<typename SeedLoc>
void search_shape(int sid, int query_block, unsigned query_iteration, char *query_buffer, char *ref_buffer, Search::Config& cfg, const HashedSeedSet* target_seeds)
{
	using SA = SeedArray<SeedLoc>;
	Partition<SeedPartition> p((SeedPartition)seedp_count(cfg.seedp_bits), cfg.index_chunks);
	log_rss();
	SequenceSet& ref_seqs = cfg.target->seqs(), &query_seqs = cfg.query->seqs();
	const SeedHistogram& ref_hst = cfg.target->hst(), query_hst = cfg.query->hst();

	for (unsigned chunk = 0; chunk < p.parts; ++chunk) {
		message_stream << "Processing query block " << query_block + 1;
		if (cfg.iterated())
			message_stream << ", query iteration " << query_iteration + 1;
		message_stream << ", reference block " << (cfg.current_ref_block + 1) << "/" << cfg.ref_blocks
			<< ", shape " << (sid + 1) << "/" << shapes.count();
		if (cfg.index_chunks > 1)
			message_stream << ", index chunk " << chunk + 1 << "/" << cfg.index_chunks;
		message_stream << '.' << endl;
		const SeedPartitionRange range(p.begin(chunk), p.end(chunk));
		current_range = range;

		TaskTimer timer("Building reference seed array", true);
		SA *ref_idx;
		const EnumCfg enum_ref{ &ref_hst.partition(), sid, sid + 1, cfg.seed_encoding, nullptr, false, false, cfg.seed_complexity_cut,
			query_seeds_bitset.get() || (bool)query_seeds_hashed ? MaskingAlgo::NONE : cfg.soft_masking,
			cfg.minimizer_window, false, false, cfg.sketch_size };
		if (query_seeds_bitset.get())
			ref_idx = new SA(*cfg.target, ref_hst.get(sid), range, cfg.seedp_bits, ref_buffer, query_seeds_bitset.get(), enum_ref);
		else if (query_seeds_hashed.get())
			ref_idx = new SA(*cfg.target, ref_hst.get(sid), range, cfg.seedp_bits, ref_buffer, query_seeds_hashed.get(), enum_ref);
			//ref_idx = new SeedArray(ref_seqs, sid, range, query_seeds_hashed.get(), true);
		else
			ref_idx = new SA(*cfg.target, ref_hst.get(sid), range, cfg.seedp_bits, ref_buffer, &no_filter, enum_ref);
		timer.finish();
		log_rss();

		timer.go("Building query seed array");
		SA* query_idx;
		EnumCfg enum_query{ target_seeds ? nullptr : &query_hst.partition(), sid, sid + 1, cfg.seed_encoding, cfg.query_skip.get(),
			false, true, cfg.seed_complexity_cut, cfg.soft_masking, cfg.minimizer_window, static_cast<bool>(query_seeds_hashed.get()),
			static_cast<bool>(query_seeds_hashed.get()), cfg.sketch_size };
		if (target_seeds)
			query_idx = new SA(*cfg.query, range, cfg.seedp_bits, target_seeds, enum_query);
		else
			query_idx = new SA(*cfg.query, query_hst.get(sid), range, cfg.seedp_bits, query_buffer, &no_filter, enum_query);
		timer.finish();
		log_rss();

		log_stream << "Indexed query seeds = " << Util::String::ratio_percentage(query_idx->size(), query_seqs.letters())
			<< ", reference seeds = " << Util::String::ratio_percentage(ref_idx->size(), ref_seqs.letters()) << endl;
		log_stream << "Soft masked letters = " << Util::String::ratio_percentage(cfg.query->soft_masked_letters(), cfg.query->seqs().letters())
			<< ", " << Util::String::ratio_percentage(cfg.target->soft_masked_letters(), cfg.target->seqs().letters()) << endl;
		/*log_stream << "Low complexity seeds = " << Util::String::ratio_percentage(query_idx->stats().low_complexity_seeds, query_idx->stats().good_seed_positions)
			<< ", " << Util::String::ratio_percentage(ref_idx->stats().low_complexity_seeds, ref_idx->stats().good_seed_positions) << endl;*/

		timer.go("Computing hash join");
		atomic<SeedPartition> seedp(0);
		vector<std::thread> threads;
		vector<DoubleArray<SeedLoc>> query_seed_hits(range.size()), ref_seed_hits(range.size());
		for (int i = 0; i < config.threads_; ++i)
			threads.emplace_back(seed_join_worker<SeedLoc>, query_idx, ref_idx, &seedp, range.size(), query_seed_hits.data(), ref_seed_hits.data());
		for (auto &t : threads)
			t.join();
		timer.finish();
		log_rss();

		if(config.freq_masking && !config.lin_stage1 && !cfg.lin_stage1_target) {
			timer.go("Building seed filter");
			frequent_seeds.build(sid, range, query_seed_hits.data(), ref_seed_hits.data(), cfg);
		}
		else
			Search::mask_seeds(shapes[sid], range, query_seed_hits.data(), ref_seed_hits.data(), cfg);

		log_rss();
		unique_ptr<KmerRanking> kmer_ranking;
		if (keep_target_id(cfg) && config.lin_stage1) {
			timer.go("Building kmer ranking");
			kmer_ranking.reset(config.kmer_ranking ? new KmerRanking(cfg.query->seqs(), range.size(), query_seed_hits.data(), ref_seed_hits.data())
				: new KmerRanking(cfg.query->seqs()));
		}

		Search::Context* context = nullptr;
		const vector<uint32_t> patterns = shapes.patterns(0, sid + 1);
		context = new Search::Context{ {patterns.data(), patterns.data() + patterns.size() - 1 },
			{patterns.data(), patterns.data() + patterns.size() },
			score_matrix.rawscore(config.short_query_ungapped_bitscore),
			kmer_ranking.get(),
			seedp_mask(cfg.seedp_bits)
		};

		timer.go("Searching alignments");
		seedp = 0;
		threads.clear();
		for (int i = 0; i < config.threads_; ++i)
			threads.emplace_back(search_worker<SeedLoc>, &seedp, range.size(), sid, i, query_seed_hits.data(), ref_seed_hits.data(), context, &cfg);
		for (auto &t : threads)
			t.join();
		timer.finish();
		log_rss();

		timer.go("Deallocating memory");
		delete ref_idx;
		delete query_idx;
		delete context;
		kmer_ranking.reset();
		timer.finish();
		log_rss();
	}
}

void search_shape(unsigned sid, int query_block, unsigned query_iteration, char* query_buffer, char* ref_buffer, Search::Config& cfg, const HashedSeedSet* target_seeds) {
	if (keep_target_id(cfg))
		search_shape<PackedLocId>(sid, query_block, query_iteration, query_buffer, ref_buffer, cfg, target_seeds);
	else
		search_shape<PackedLoc>(sid, query_block, query_iteration, query_buffer, ref_buffer, cfg, target_seeds);
}

}
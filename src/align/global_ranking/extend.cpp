/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

#include <unordered_map>
#include <mutex>
#include <memory>
#include <thread>
#include <atomic>
#include "global_ranking.h"
#include "output/output.h"
#include "align/target.h"
#include "data/queries.h"
#include "masking/masking.h"
#include "data/sequence_file.h"
#include "util/parallel/simple_thread_pool.h"

using std::unique_ptr;
using std::mutex;
using std::thread;
using std::unordered_map;
using std::endl;
using std::pair;
using std::vector;

namespace Extension { namespace GlobalRanking {

typedef unordered_map<OId, BlockId> TargetMap;

void extend_query(const QueryList& query_list, const TargetMap& db2block_id, const Search::Config& cfg, Statistics& stats, std::pmr::monotonic_buffer_resource& pool) {
	SeedHitList l;
	const size_t n = query_list.targets.size();
	l.target_block_ids.reserve(n);
	l.target_scores.reserve(n);
	l.seed_hits.reserve(n, 0);
	for (size_t i = 0; i < n; ++i) {
		l.target_block_ids.push_back(db2block_id.at(query_list.targets[i].database_id));
		l.target_scores.push_back({ (uint32_t)i, query_list.targets[i].score });
		l.seed_hits.next();
		l.seed_hits.push_back({ 0,0,query_list.targets[i].score,0 });
	}
	
	DP::Flags flags = DP::Flags::FULL_MATRIX;
	
	vector<Match> matches = Extension::extend(
		query_list.query_block_id,
		cfg,
		stats,
		flags,
		l,
		pool);

	TextBuffer* buf = Extension::generate_output(matches, query_list.query_block_id, stats, cfg);
	if (!matches.empty() && (!config.unaligned.empty() || !config.aligned_file.empty())) {
		std::lock_guard<std::mutex> lock(query_aligned_mtx);
		query_aligned[query_list.query_block_id] = true;
	}
	output_sink->push(query_list.query_block_id, buf);
}

void align_worker(const std::atomic<bool>& stop, File* query_list, const TargetMap* db2block_id, const Search::Config* cfg, uint32_t* next_query) {
	std::pmr::monotonic_buffer_resource pool;
	QueryList input;
	Statistics stats;
	while (!stop.load(std::memory_order_relaxed) && (input = fetch_query_targets(*query_list, *next_query), !input.targets.empty())) {
		for (uint32_t i = input.last_query_block_id; i < input.query_block_id; ++i)
			output_sink->push(i, nullptr);
		extend_query(input, *db2block_id, *cfg, stats, pool);
	}
	statistics += stats;
}

void extend(SequenceFile& db, File& merged_query_list, BitVector& ranking_db_filter, Search::Config& cfg, File& master_out) {
	TaskTimer timer("Loading reference sequences");
	db.set_seqinfo_ptr(0);
	db.flags() |= SequenceFile::Flags::SEQS;
	cfg.target.reset(db.load_seqs(INT64_MAX, 0, &ranking_db_filter));
	TargetMap db2block_id;
	const BlockId db_count = cfg.target->seqs().size();
	db2block_id.reserve(db_count);
	for (BlockId i = 0; i < db_count; ++i)
		db2block_id[cfg.target->block_id2oid(i)] = i;
	timer.finish();
	*log_stream << "#Ranked database sequences: " << cfg.target->seqs().size() << endl;

	if (cfg.target_masking != MaskingAlgo::NONE) {
		timer.go("Masking reference");
		const MaskingStat stats = mask_seqs(cfg.target->seqs(), Masking::get(), true, cfg.target_masking);
		timer.finish();
		stats.print(*log_stream);
	}

	timer.go("Computing alignments");
	OutputWriter writer{ &master_out };
	output_sink.reset(new ReorderQueue<TextBuffer*, OutputWriter>(0, writer));
	uint32_t next_query = 0;
	SimpleThreadPool pool;
	for (size_t i = 0; i < (config.threads_align ? config.threads_align : config.threads_); ++i)
		pool.spawn(align_worker, &merged_query_list, &db2block_id, &cfg, &next_query);
	pool.join_all();

	timer.go("Cleaning up");
	merged_query_list.close();
	output_sink.reset();
	cfg.target.reset();
}

void extend_query(BlockId source_query_block_id, const TargetMap& db2block_id, Search::Config& cfg, Statistics& stats, std::pmr::monotonic_buffer_resource& pool) {
	const size_t N = config.global_ranking_targets;
	SeedHitList l;
	vector<Hit>::const_iterator table_begin = cfg.ranking_table->cbegin() + source_query_block_id * N, table_end = table_begin + N;
	while (table_end > table_begin && (table_end - 1)->score == 0) --table_end;
	const size_t n = table_end - table_begin;
	TextBuffer* buf = nullptr;
	if (n) {
		l.target_block_ids.reserve(n);
		l.target_scores.reserve(n);
		l.seed_hits.reserve(n, 0);
		for (size_t i = 0; i < n; ++i) {
			l.target_block_ids.push_back(db2block_id.at(table_begin[i].oid));
			l.target_scores.push_back({ (uint32_t)i, table_begin[i].score });
			l.seed_hits.next();
			l.seed_hits.push_back({ 0,0,table_begin[i].score, table_begin[i].context });
		}

		DP::Flags flags = DP::Flags::FULL_MATRIX;

		vector<Match> matches = Extension::extend(
			source_query_block_id,
			cfg,
			stats,
			flags,
			l,
			pool);

		buf = cfg.iterated() ? Extension::generate_intermediate_output(matches, source_query_block_id, cfg) : Extension::generate_output(matches, source_query_block_id, stats, cfg);

		if (!matches.empty() && cfg.track_aligned_queries) {
			std::lock_guard<std::mutex> lock(query_aligned_mtx);
			if (!query_aligned[source_query_block_id]) {
				query_aligned[source_query_block_id] = true;
				++cfg.iteration_query_aligned;
			}
		}
	}
	output_sink->push(source_query_block_id, buf);
}


static BitVector db_filter(const Search::Config::RankingTable& table, size_t db_size) {
	BitVector v(db_size);
	for (const Hit& hit : table)
		if (hit.score)
			v.set(hit.oid);
	return v;
}

void extend(Search::Config& cfg, File& out) {
	TaskTimer timer("Listing target sequences");
		const BitVector filter = db_filter(*cfg.ranking_table, cfg.db->sequence_count().value());
	timer.go("Loading target sequences");
	cfg.db->set_seqinfo_ptr(0);
	cfg.db->flags() |= SequenceFile::Flags::SEQS;
	if (!flag_any(cfg.db->format_flags(), SequenceFile::FormatFlags::TITLES_LAZY))
		cfg.db->flags() |= SequenceFile::Flags::TITLES;
	if (bool(cfg.output_format->flags & Output::Flags::FULL_TITLES))
		cfg.db->flags() |= SequenceFile::Flags::FULL_TITLES;
	cfg.target.reset(cfg.db->load_seqs(INT64_MAX, 0, &filter));
	TargetMap db2block_id;
	const BlockId db_count = cfg.target->seqs().size();
	db2block_id.reserve(db_count);
	for (BlockId i = 0; i < db_count; ++i)
		db2block_id[cfg.target->block_id2oid(i)] = i;
	timer.finish();
	*log_stream << "#Ranked database sequences: " << db_count << endl;

	if (cfg.target_masking != MaskingAlgo::NONE) {
		timer.go("Masking reference");
		const MaskingStat stats = mask_seqs(cfg.target->seqs(), Masking::get(), true, cfg.target_masking);
		timer.finish();
		stats.print(*log_stream);
	}

	if (cfg.iterated()) {
		cfg.current_ref_block = 0;
		cfg.db->init_dict_block(0, db_count, true);
	}
	else
		cfg.db->init_random_access(cfg.current_query_block, 0);

	timer.go("Computing alignments");
	OutputWriter writer{ &out };
	output_sink.reset(new ReorderQueue<TextBuffer*, OutputWriter>(0, writer));

	std::atomic<BlockId> next_query(0);
	const BlockId query_count = cfg.query->seqs().size() / align_mode.query_contexts;
	auto worker = [&next_query, &db2block_id, &cfg, query_count](const std::atomic<bool>& stop) {
		Statistics stats;
		BlockId q;
		std::pmr::monotonic_buffer_resource pool;
		while (!stop.load(std::memory_order_relaxed) && (q = next_query++) < query_count) extend_query(q, db2block_id, cfg, stats, pool);
		statistics += stats;
		};

	SimpleThreadPool pool;
	for (size_t i = 0; i < (config.threads_align ? config.threads_align : config.threads_); ++i)
		pool.spawn(worker);
	pool.join_all();

	timer.go("Deallocating memory");
	cfg.target.reset();
	output_sink.reset();
	if (!cfg.iterated())
		cfg.db->end_random_access();
	else {
		cfg.db->close_dict_block(true);
		IntermediateRecord::finish_file(out);
	}
}

}}
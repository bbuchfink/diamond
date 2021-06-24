/****
DIAMOND protein aligner
Copyright (C) 2020 Max Planck Society for the Advancement of Science e.V.

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

#include <unordered_map>
#include <mutex>
#include <memory>
#include <thread>
#include "global_ranking.h"
#include "../output/output.h"
#include "../target.h"
#include "../dp/dp.h"
#include "../data/queries.h"
#include "../basic/masking.h"

using std::unique_ptr;
using std::mutex;
using std::thread;
using std::unordered_map;
using std::endl;

namespace Extension { namespace GlobalRanking {

typedef unordered_map<uint32_t, uint32_t> TargetMap;

void extend_query(const QueryList& query_list, const TargetMap& db2block_id, const Search::Config& cfg, Statistics& stats) {
	thread_local vector<uint32_t> target_block_ids;
	thread_local vector<TargetScore> target_scores;
	thread_local FlatArray<SeedHit> seed_hits;
	const size_t n = query_list.targets.size();
	target_block_ids.clear();
	target_block_ids.reserve(n);
	target_scores.clear();
	target_scores.reserve(n);
	seed_hits.clear();
	seed_hits.reserve(n);
	for (size_t i = 0; i < n; ++i) {
		target_block_ids.push_back(db2block_id.at(query_list.targets[i].database_id));
		target_scores.push_back({ (uint32_t)i, query_list.targets[i].score });
		seed_hits.next();
		seed_hits.push_back({ 0,0,query_list.targets[i].score,0 });
	}
	
	int flags = DP::FULL_MATRIX;
	
	vector<Match> matches = Extension::extend(
		query_list.query_block_id,
		cfg,
		stats,
		flags,
		seed_hits,
		target_block_ids,
		target_scores);

	TextBuffer* buf = Extension::generate_output(matches, query_list.query_block_id, stats, cfg);
	if (!matches.empty() && (!config.unaligned.empty() || !config.aligned_file.empty())) {
		std::lock_guard<std::mutex> lock(query_aligned_mtx);
		query_aligned[query_list.query_block_id] = true;
	}
	OutputSink::get().push(query_list.query_block_id, buf);
}

void align_worker(InputFile* query_list, const TargetMap* db2block_id, const Search::Config* cfg, uint32_t* next_query) {
	try {
		QueryList input;
		Statistics stats;
		while (input = fetch_query_targets(*query_list, *next_query), !input.targets.empty()) {
			for (uint32_t i = input.last_query_block_id; i < input.query_block_id; ++i)
				OutputSink::get().push(i, nullptr);
			extend_query(input, *db2block_id, *cfg, stats);
		}
		statistics += stats;
	}
	catch (std::exception& e) {
		exit_with_error(e);
	}
}

void extend(SequenceFile& db, TempFile& merged_query_list, BitVector& ranking_db_filter, Search::Config& cfg, Consumer& master_out) {
	task_timer timer("Loading reference sequences");
	InputFile query_list(merged_query_list);
	db.set_seqinfo_ptr(0);
	cfg.target.reset(db.load_seqs(SIZE_MAX, false, &ranking_db_filter, true));
	TargetMap db2block_id;
	const size_t db_count = cfg.target->seqs().size();
	db2block_id.reserve(db_count);
	for (size_t i = 0; i < db_count; ++i)
		db2block_id[cfg.target->block_id2oid(i)] = i;
	timer.finish();
	verbose_stream << "#Ranked database sequences: " << cfg.target->seqs().size() << endl;

	if (config.masking == 1) {
		timer.go("Masking reference");
		size_t n = mask_seqs(cfg.target->seqs(), Masking::get());
		timer.finish();
		log_stream << "Masked letters: " << n << endl;
	}

	timer.go("Computing alignments");
	OutputSink::instance.reset(new OutputSink(0, &master_out));
	uint32_t next_query = 0;
	vector<thread> threads;
	for (size_t i = 0; i < (config.threads_align ? config.threads_align : config.threads_); ++i)
		threads.emplace_back(align_worker, &query_list, &db2block_id, &cfg, &next_query);
	for (auto& i : threads)
		i.join();

	timer.go("Cleaning up");
	query_list.close_and_delete();
	cfg.target.reset();
}

void extend_query(size_t source_query_block_id, const TargetMap& db2block_id, Search::Config& cfg, Statistics& stats) {
	const size_t N = config.global_ranking_targets;
	thread_local vector<uint32_t> target_block_ids;
	thread_local vector<TargetScore> target_scores;
	thread_local FlatArray<SeedHit> seed_hits;
	vector<Hit>::const_iterator table_begin = cfg.ranking_table->cbegin() + source_query_block_id * N, table_end = table_begin + N;
	while (table_end > table_begin && (table_end - 1)->score == 0) --table_end;
	const size_t n = table_end - table_begin;
	TextBuffer* buf = nullptr;
	if (n) {
		target_block_ids.clear();
		target_block_ids.reserve(n);
		target_scores.clear();
		target_scores.reserve(n);
		seed_hits.clear();
		seed_hits.reserve(n);
		for (size_t i = 0; i < n; ++i) {
			target_block_ids.push_back(db2block_id.at(table_begin[i].oid));
			target_scores.push_back({ (uint32_t)i, table_begin[i].score });
			seed_hits.next();
			seed_hits.push_back({ 0,0,table_begin[i].score, table_begin[i].context });
		}

		int flags = DP::FULL_MATRIX;

		vector<Match> matches = Extension::extend(
			source_query_block_id,
			cfg,
			stats,
			flags,
			seed_hits,
			target_block_ids,
			target_scores);

		buf = cfg.iterated() ? Extension::generate_intermediate_output(matches, source_query_block_id, cfg) : Extension::generate_output(matches, source_query_block_id, stats, cfg);

		if (!matches.empty() && cfg.track_aligned_queries) {
			std::lock_guard<std::mutex> lock(query_aligned_mtx);
			if (!query_aligned[source_query_block_id]) {
				query_aligned[source_query_block_id] = true;
				++cfg.iteration_query_aligned;
			}
		}
	}
	OutputSink::get().push(source_query_block_id, buf);
}


static BitVector db_filter(const Search::Config::RankingTable& table, size_t db_size) {
	BitVector v(db_size);
	for (const Hit& hit : table)
		if (hit.score)
			v.set(hit.oid);
	return v;
}

void extend(Search::Config& cfg, Consumer& out) {
	task_timer timer("Listing target sequences");
	const BitVector filter = db_filter(*cfg.ranking_table, cfg.db->sequence_count());
	timer.go("Loading target sequences");
	cfg.db->set_seqinfo_ptr(0);
	cfg.target.reset(cfg.db->load_seqs(SIZE_MAX, cfg.db->load_titles() == SequenceFile::LoadTitles::SINGLE_PASS, &filter, true));
	TargetMap db2block_id;
	const size_t db_count = cfg.target->seqs().size();
	db2block_id.reserve(db_count);
	for (size_t i = 0; i < db_count; ++i)
		db2block_id[cfg.target->block_id2oid(i)] = i;
	timer.finish();
	verbose_stream << "#Ranked database sequences: " << db_count << endl;

	if (config.masking == 1 || config.target_seg == 1) {
		timer.go("Masking reference");
		size_t n;
		if (config.target_seg == 1)
			n = mask_seqs(cfg.target->seqs(), Masking::get(), true, Masking::Algo::SEG);
		else
			n = mask_seqs(cfg.target->seqs(), Masking::get());
		timer.finish();
		log_stream << "Masked letters: " << n << endl;
	}

	if (cfg.iterated()) {
		current_ref_block = 0;
		cfg.db->init_dict_block(0, db_count, true);
	}
	else
		cfg.db->init_random_access(current_query_chunk, 0);

	timer.go("Computing alignments");
	OutputSink::instance.reset(new OutputSink(0, &out));
	std::atomic_size_t next_query(0);
	const size_t query_count = cfg.query->seqs().size() / align_mode.query_contexts;
	auto worker = [&next_query, &db2block_id, &cfg, query_count] {
		try {
			Statistics stats;
			size_t q;
			while ((q = next_query++) < query_count) extend_query(q, db2block_id, cfg, stats);
			statistics += stats;
		}
		catch (std::exception& e) {
			exit_with_error(e);
		}
	};

	vector<thread> threads;
	for (size_t i = 0; i < (config.threads_align ? config.threads_align : config.threads_); ++i)
		threads.emplace_back(worker);
	for (auto& i : threads)
		i.join();

	timer.go("Deallocating memory");
	cfg.target.reset();
	if (!cfg.iterated())
		cfg.db->end_random_access();
	else {
		cfg.db->close_dict_block(true);
		IntermediateRecord::finish_file(out);
	}
}

}}
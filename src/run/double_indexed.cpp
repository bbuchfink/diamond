/****
DIAMOND protein aligner
Copyright (C) 2013-2021 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen

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

#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <algorithm>
#include <cstdio>
#include "../data/reference.h"
#include "../data/queries.h"
#include "../basic/statistics.h"
#include "../basic/shape_config.h"
#include "../util/seq_file_format.h"
#include "../data/load_seqs.h"
#include "../output/output_format.h"
#include "../data/frequent_seeds.h"
#include "../output/daa/daa_write.h"
#include "../data/taxonomy.h"
#include "../basic/masking.h"
#include "../data/block.h"
#include "../search/search.h"
#include "workflow.h"
#include "../util/io/consumer.h"
#include "../util/parallel/thread_pool.h"
#include "../util/parallel/multiprocessing.h"
#include "../util/parallel/parallelizer.h"
#include "../util/system/system.h"
#include "../align/target.h"
#include "../data/seed_set.h"
#include "../util/data_structures/deque.h"
#include "../align/global_ranking/global_ranking.h"
#include "../align/align.h"
#include "config.h"

using std::unique_ptr;
using std::endl;
using std::list;
using std::shared_ptr;

namespace Search {

static const size_t MAX_INDEX_QUERY_SIZE = 32 * MEGABYTES;
static const size_t MAX_HASH_SET_SIZE = 8 * MEGABYTES;
static const size_t MIN_QUERY_INDEXED_DB_SIZE = 256 * MEGABYTES;

static const string label_align = "align";
static const string stack_align_todo = label_align + "_todo";
static const string stack_align_wip = label_align + "_wip";
static const string stack_align_done = label_align + "_done";

static const string label_join = "join";
static const string stack_join_todo = label_join + "_todo";
static const string stack_join_wip = label_join + "_wip";
static const string stack_join_redo = label_join + "_redo";
static const string stack_join_done = label_join + "_done";

static bool use_query_index(const size_t table_size) {
	return table_size <= std::max(MAX_HASH_SET_SIZE, l3_cache_size());
}

string get_ref_part_file_name(const string & prefix, size_t query, string suffix="") {
	if (suffix.size() > 0)
		suffix.append("_");
	const string file_name = append_label(prefix + "_" + suffix, query);
	return join_path(config.parallel_tmpdir, file_name);
}

string get_ref_block_tmpfile_name(size_t query, size_t block) {
	const string file_name = append_label("ref_block_", query) + append_label("_", block);
	return join_path(config.parallel_tmpdir, file_name);
}

void run_ref_chunk(SequenceFile &db_file,
	const unsigned query_chunk,
	const unsigned query_iteration,
	char *query_buffer,
	Consumer &master_out,
	PtrVector<TempFile> &tmp_file,
	Config& cfg)
{
	log_rss();
	auto& ref_seqs = cfg.target->seqs();
	auto& query_seqs = cfg.query->seqs();

	if (config.comp_based_stats == Stats::CBS::COMP_BASED_STATS_AND_MATRIX_ADJUST || flag_any(output_format->flags, Output::Flags::TARGET_SEQS))
		cfg.target->unmasked_seqs() = ref_seqs;

	task_timer timer;
	if (config.masking == 1 && !config.no_ref_masking && !cfg.lazy_masking) {
		timer.go("Masking reference");
		size_t n = mask_seqs(ref_seqs, Masking::get());
		timer.finish();
		log_stream << "Masked letters: " << n << endl;
	}

	const bool daa = *output_format == Output_format::daa;
	const bool persist_dict = daa || cfg.iterated();
	if(((blocked_processing || daa) && !config.global_ranking_targets) || cfg.iterated()) {
		timer.go("Initializing dictionary");
		if (config.multiprocessing || (current_ref_block == 0 && (!daa || query_chunk == 0) && query_iteration == 0))
			db_file.init_dict(query_chunk, current_ref_block);
		if(!config.global_ranking_targets)
			db_file.init_dict_block(current_ref_block, ref_seqs.size(), persist_dict);
	}

	timer.go("Initializing temporary storage");
	if (config.global_ranking_targets)
		;// cfg.global_ranking_buffer.reset(new Config::RankingBuffer());
	else
		cfg.seed_hit_buf.reset(new AsyncBuffer<Search::Hit>(query_seqs.size() / align_mode.query_contexts,
			config.tmpdir,
			cfg.query_bins,
			{ cfg.target->long_offsets(), align_mode.query_contexts }));

	if (!config.swipe_all) {
		timer.go("Building reference histograms");
		if(query_seeds_bitset.get())
			cfg.target->hst() = Partitioned_histogram(ref_seqs, true, query_seeds_bitset.get(), cfg.seed_encoding, nullptr);
		else if (query_seeds_hashed.get())
			cfg.target->hst() = Partitioned_histogram(ref_seqs, true, query_seeds_hashed.get(), cfg.seed_encoding, nullptr);
		else
			cfg.target->hst() = Partitioned_histogram(ref_seqs, false, &no_filter, cfg.seed_encoding, nullptr);

		timer.go("Allocating buffers");
		char *ref_buffer = SeedArray::alloc_buffer(cfg.target->hst(), cfg.index_chunks);
		timer.finish();

		HashedSeedSet* target_seeds = nullptr;
		if (config.target_indexed) {
			timer.go("Loading database seed index");
			target_seeds = new HashedSeedSet(db_file.file_name() + ".seed_idx");
			timer.finish();
		}

		for (unsigned i = 0; i < shapes.count(); ++i) {
			if(config.global_ranking_targets)
				cfg.global_ranking_buffer.reset(new Config::RankingBuffer());
			search_shape(i, query_chunk, query_iteration, query_buffer, ref_buffer, cfg, target_seeds);
			if (config.global_ranking_targets)
				Extension::GlobalRanking::update_table(cfg);
		}

		timer.go("Deallocating buffers");
		delete[] ref_buffer;
		delete target_seeds;

		timer.go("Clearing query masking");
		Frequent_seeds::clear_masking(query_seqs);
	}

	Consumer* out;
	const bool temp_output = (blocked_processing || cfg.iterated()) && !config.global_ranking_targets;
	if (temp_output) {
		timer.go("Opening temporary output file");
		if (config.multiprocessing) {
			const string file_name = get_ref_block_tmpfile_name(query_chunk, current_ref_block);
			tmp_file.push_back(new TempFile(file_name));
		} else {
			tmp_file.push_back(new TempFile());
		}
		out = &tmp_file.back();
	}
	else
		out = &master_out;

	if (config.target_seg == 1 && !cfg.lazy_masking) {
		timer.go("SEG masking targets");
		mask_seqs(ref_seqs, Masking::get(), true, Masking::Algo::SEG);
	}

	if (config.global_ranking_targets) {
		/*timer.go("Updating ranking table");
		Extension::GlobalRanking::update_table(cfg);
		cfg.global_ranking_buffer.reset();*/
	}
	else {
		timer.go("Computing alignments");
		align_queries(out, cfg);
		cfg.seed_hit_buf.reset();
	}	

	if (temp_output)
		IntermediateRecord::finish_file(*out);

	timer.go("Deallocating reference");
	cfg.target.reset();
	cfg.db->close_dict_block(persist_dict);
	timer.finish();
}

void run_query_iteration(const unsigned query_chunk,
	const unsigned query_iteration,
	Consumer& master_out,
	OutputFile* unaligned_file,
	OutputFile* aligned_file,
	PtrVector<TempFile>& tmp_file,
	Config& options)
{
	task_timer timer;
	auto P = Parallelizer::get();
	auto& query_seqs = options.query->seqs();
	auto& query_ids = options.query->ids();
	auto& db_file = *options.db;
	if (query_iteration > 0)
		options.query_skip.reset(new vector<bool> (query_aligned));

	if (config.algo == ::Config::Algo::AUTO &&
		(!sensitivity_traits.at(config.sensitivity).support_query_indexed
			|| query_seqs.letters() > MAX_INDEX_QUERY_SIZE
			|| options.db->letters() < MIN_QUERY_INDEXED_DB_SIZE
			|| config.target_indexed))
		config.algo = ::Config::Algo::DOUBLE_INDEXED;
	if (config.algo == ::Config::Algo::AUTO || config.algo == ::Config::Algo::QUERY_INDEXED) {
		timer.go("Building query seed set");
		query_seeds_hashed.reset(new HashedSeedSet(query_seqs, options.query_skip.get()));
		if (config.algo == ::Config::Algo::AUTO && !use_query_index(query_seeds_hashed->max_table_size())) {
			config.algo = ::Config::Algo::DOUBLE_INDEXED;
			query_seeds_hashed.reset();
		}
		else {
			config.algo = ::Config::Algo::QUERY_INDEXED;
			options.seed_encoding = SeedEncoding::HASHED;
		}
		timer.finish();
	}
	if (config.algo == ::Config::Algo::CTG_SEED) {
		timer.go("Building query seed set");
		query_seeds_bitset.reset(new SeedSet(query_seqs, 1.0, options.query_skip.get()));
		options.seed_encoding = SeedEncoding::CONTIGUOUS;
		timer.finish();
	}

	::Config::set_option(options.index_chunks, config.lowmem_, 0u, config.algo == ::Config::Algo::DOUBLE_INDEXED ? sensitivity_traits.at(options.sensitivity[query_iteration]).index_chunks : 1u);
	options.lazy_masking = config.algo != ::Config::Algo::DOUBLE_INDEXED && (config.masking == 1 || config.target_seg == 1) && config.frame_shift == 0;

	options.cutoff_gapped1 = { config.gapped_filter_evalue1 };
	options.cutoff_gapped2 = { options.gapped_filter_evalue };
	options.cutoff_gapped1_new = { config.gapped_filter_evalue1 };
	options.cutoff_gapped2_new = { options.gapped_filter_evalue };

	if (current_query_chunk == 0 && query_iteration == 0) {
		message_stream << "Algorithm: " << to_string(config.algo) << endl;
		verbose_stream << "Seed frequency SD: " << options.freq_sd << endl;
		verbose_stream << "Shape configuration: " << ::shapes << endl;
	}	

	if (config.global_ranking_targets) {
		timer.go("Allocating global ranking table");
		options.ranking_table.reset(new Search::Config::RankingTable(query_seqs.size() * config.global_ranking_targets / align_mode.query_contexts));
	}

	char* query_buffer = nullptr;
	if (!config.swipe_all && !config.target_indexed) {
		timer.go("Building query histograms");
		options.query->hst() = Partitioned_histogram(query_seqs, false, &no_filter, options.seed_encoding, options.query_skip.get());

		timer.go("Allocating buffers");
		query_buffer = SeedArray::alloc_buffer(options.query->hst(), options.index_chunks);
		timer.finish();
	}

	log_rss();
	db_file.set_seqinfo_ptr(0);
	bool mp_last_chunk = false;

	if (config.multiprocessing) {
		P->create_stack_from_file(stack_align_todo, get_ref_part_file_name(stack_align_todo, query_chunk));
		auto work = P->get_stack(stack_align_todo);
		P->create_stack_from_file(stack_align_wip, get_ref_part_file_name(stack_align_wip, query_chunk));
		auto wip = P->get_stack(stack_align_wip);
		P->create_stack_from_file(stack_align_done, get_ref_part_file_name(stack_align_done, query_chunk));
		auto done = P->get_stack(stack_align_done);
		P->create_stack_from_file(stack_join_todo, get_ref_part_file_name(stack_join_todo, query_chunk));
		auto join_work = P->get_stack(stack_join_todo);

		string buf;

		while ((!file_exists("stop")) && (work->pop(buf))) {
			wip->push(buf);

			Chunk chunk = to_chunk(buf);

			P->log("SEARCH BEGIN " + std::to_string(query_chunk) + " " + std::to_string(chunk.i));

			options.target.reset(db_file.load_seqs((size_t)(0), db_file.load_titles() == SequenceFile::LoadTitles::SINGLE_PASS, options.db_filter.get(), true, options.lazy_masking, chunk));
			run_ref_chunk(db_file, query_chunk, query_iteration, query_buffer, master_out, tmp_file, options);

			tmp_file.back().close();

			size_t size_after_push = 0;
			done->push(buf, size_after_push);
			if (size_after_push == db_file.get_n_partition_chunks()) {
				join_work->push("TOKEN");
			}
			wip->remove(buf);

			P->log("SEARCH END " + std::to_string(query_chunk) + " " + std::to_string(chunk.i));
			log_rss();
		}

		tmp_file.clear();
		P->delete_stack(stack_align_todo);
		P->delete_stack(stack_align_wip);
		P->delete_stack(stack_align_done);
	}
	else {
		for (current_ref_block = 0; ; ++current_ref_block) {
			options.target.reset(db_file.load_seqs((size_t)(config.chunk_size * 1e9), db_file.load_titles() == SequenceFile::LoadTitles::SINGLE_PASS, options.db_filter.get(), true, options.lazy_masking));
			if (options.target->empty()) break;
			run_ref_chunk(db_file, query_chunk, query_iteration, query_buffer, master_out, tmp_file, options);
		}
		log_rss();
	}

	timer.go("Deallocating buffers");
	delete[] query_buffer;
	query_seeds_hashed.reset();
	query_seeds_bitset.reset();
	options.query_skip.reset();
	delete Extension::memory;

	if (config.global_ranking_targets) {
		timer.go("Computing alignments");
		Consumer* out;
		if (options.iterated()) {
			tmp_file.push_back(new TempFile());
			out = &tmp_file.back();
		}
		else {
			out = &master_out;
		}
		Extension::GlobalRanking::extend(options, *out);
		options.ranking_table.reset();
	}
}

void run_query_chunk(const unsigned query_chunk,
	Consumer &master_out,
	OutputFile *unaligned_file,
	OutputFile *aligned_file,
	Config &options)
{
	auto P = Parallelizer::get();
	task_timer timer;
	auto& db_file = *options.db;
	auto& query_seqs = options.query->seqs();
	auto& query_ids = options.query->ids();
	
	options.db_seqs = db_file.sequence_count();
	options.db_letters = db_file.letters();
	options.ref_blocks = db_file.total_blocks();

	PtrVector<TempFile> tmp_file;
	if (options.track_aligned_queries) {
		query_aligned.clear();
		query_aligned.insert(query_aligned.end(), query_ids.size(), false);
	}
	if(config.query_memory)
		Extension::memory = new Extension::Memory(query_ids.size());

	log_rss();

	size_t aligned = 0;
	for (unsigned query_iteration = 0; query_iteration < options.sensitivity.size() && aligned < query_ids.size(); ++query_iteration) {
		setup_search(options.sensitivity[query_iteration], options);
		run_query_iteration(query_chunk, query_iteration, master_out, unaligned_file, aligned_file, tmp_file, options);
		if (options.iterated()) {
			aligned += options.iteration_query_aligned;
			message_stream << "Aligned " << options.iteration_query_aligned << '/' << query_ids.size() << " queries in this iteration, "
				<< aligned << '/' << query_ids.size() << " total." << endl;
			options.iteration_query_aligned = 0;
		}
	}

	log_rss();

	if (!config.single_chunk && (blocked_processing || config.multiprocessing || options.iterated())) {
		if(!config.global_ranking_targets) timer.go("Joining output blocks");

		if (config.multiprocessing) {
			P->create_stack_from_file(stack_join_todo, get_ref_part_file_name(stack_join_todo, query_chunk));
			auto work = P->get_stack(stack_join_todo);

			string buf;

			if ((! file_exists("stop")) && (work->pop(buf) > 0)) {
				P->log("JOIN BEGIN "+std::to_string(query_chunk));

				P->create_stack_from_file(stack_join_wip, get_ref_part_file_name(stack_join_wip, query_chunk));
				auto wip = P->get_stack(stack_join_wip);
				wip->clear();
				P->create_stack_from_file(stack_join_done, get_ref_part_file_name(stack_join_done, query_chunk));
				auto done = P->get_stack(stack_join_done);
				done->clear();

				wip->push(buf);
				work->clear();

				current_ref_block = config.join_chunks > 0 ? config.join_chunks : db_file.get_n_partition_chunks();

				vector<string> tmp_file_names;
				for (size_t i=0; i<current_ref_block; ++i) {
					tmp_file_names.push_back(get_ref_block_tmpfile_name(query_chunk, i));
				}

				const string query_chunk_output_file = append_label(config.output_file + "_", current_query_chunk);
				Consumer *query_chunk_out(new OutputFile(query_chunk_output_file, config.compressor()));
				// if (*output_format != Output_format::daa)
				// 	output_format->print_header(*query_chunk_out, align_mode.mode, config.matrix.c_str(), score_matrix.gap_open(), score_matrix.gap_extend(), config.max_evalue, query_ids::get()[0],
				// 		unsigned(align_mode.query_translated ? query_source_seqs::get()[0].length() : query_seqs::get()[0].length()));

				join_blocks(current_ref_block, *query_chunk_out, tmp_file, options, db_file, tmp_file_names);

				// if (*output_format == Output_format::daa)
				// 	// finish_daa(*static_cast<OutputFile*>(query_chunk_out), *db_file);
				// 	throw std::runtime_error("output_format::daa");
				// else
				// 	output_format->print_footer(*query_chunk_out);
				query_chunk_out->finalize();
				delete query_chunk_out;

				done->push(buf);
				wip->pop(buf);

				for (auto f : tmp_file_names) {
					std::remove(f.c_str());
				}

				P->delete_stack(stack_join_wip);
				P->delete_stack(stack_join_done);

				P->log("JOIN END "+std::to_string(query_chunk));
			}
			P->delete_stack(stack_join_todo);
		} else {
			if (!tmp_file.empty())
				join_blocks(current_ref_block, master_out, tmp_file, options, db_file);
		}
	}

	if (unaligned_file) {
		timer.go("Writing unaligned queries");
		write_unaligned(*options.query, unaligned_file);
	}
	if (aligned_file) {
		timer.go("Writing aligned queries");
		write_aligned(*options.query, aligned_file);
	}

	timer.go("Deallocating queries");
	options.query.reset();
}

void master_thread(task_timer &total_timer, Config &options)
{
	SequenceFile* db_file = options.db.get();
	if (config.multiprocessing && config.mp_recover) {
		const size_t max_assumed_query_chunks = 65536;
		for (size_t i=0; i<max_assumed_query_chunks; ++i) {
			const string file_align_todo = get_ref_part_file_name(stack_align_todo, i);
			if (! file_exists(file_align_todo)) {
				break;
			} else {
				auto stack_todo = FileStack(file_align_todo);
				const string file_align_wip = get_ref_part_file_name(stack_align_wip, i);
				auto stack_wip = FileStack(file_align_wip);
				string buf;
				int j = 0;
				while (stack_wip.pop_non_locked(buf)) {
					stack_todo.push_non_locked(buf);
					j++;
				}
				if (j > 0)
					message_stream << "Restored " << j << " align chunks for query " << i << endl;
			}
			const string file_join_wip = get_ref_part_file_name(stack_join_wip, i);
			auto stack_wip = FileStack(file_join_wip);
			if (stack_wip.size() > 0) {
				const string file_join_todo = get_ref_part_file_name(stack_join_todo, i);
				auto stack_todo = FileStack(file_join_todo);
				string buf;
				int j = 0;
				while (stack_wip.pop_non_locked(buf)) {
					stack_todo.push_non_locked(buf);
					j++;
				}
				if (j > 0)
					message_stream << "Restored join of query " << i << endl;
			}
		}
		if (file_exists("stop")) {
			std::remove("stop");
			message_stream << "Removed \'stop\' file" << endl;
		}
		return;
	}

	auto P = Parallelizer::get();
	if (config.multiprocessing) {
		P->init(config.parallel_tmpdir);
		if (config.join_chunks == 0) {
			db_file->create_partition_balanced((size_t)(config.chunk_size*1e9));
		}
	}

	task_timer timer("Opening the input file", true);
	const Sequence_file_format *format_n = nullptr;
	bool paired_mode = false;
	if (!options.self) {
		if (config.query_file.empty() && !options.query_file)
			std::cerr << "Query file parameter (--query/-q) is missing. Input will be read from stdin." << endl;
		if (!options.query_file) {
			options.query_file.reset(new list<TextInputFile>);
			for(const string& f : config.query_file)
				options.query_file->emplace_back(f);
			if (options.query_file->empty())
				options.query_file->emplace_back("");
			paired_mode = options.query_file->size() == 2;
		}
		format_n = guess_format(options.query_file->front());
	}

	current_query_chunk = 0;
	size_t query_file_offset = 0;

	if (config.multiprocessing && config.mp_init) {
		task_timer timer("Counting query blocks", true);

		size_t block_count = 0;
		do {
			if (options.self) {
				db_file->set_seqinfo_ptr(query_file_offset);
				options.query.reset(db_file->load_seqs((size_t)(config.chunk_size * 1e9), true, options.db_filter.get()));
				query_file_offset = db_file->tell_seq();
			} else {
				options.query.reset(new Block(options.query_file->begin(), options.query_file->end(), *format_n, (size_t)(config.chunk_size * 1e9), input_value_traits, config.store_query_quality, false, paired_mode ? 2 : 1));
			}
			++block_count;
		} while (!options.query->empty());
		if (options.self) {
			db_file->set_seqinfo_ptr(0);
			query_file_offset = 0;
		} else {
			options.query_file->front().rewind();
		}

		for (size_t i = 0; i < block_count - 1; ++i) {
			const string annotation = "# query_chunk=" + std::to_string(i);
			db_file->save_partition(get_ref_part_file_name(stack_align_todo, i), annotation);
		}

		return;
	}

	current_query_chunk = 0;

	timer.go("Opening the output file");
	if (!options.out)
		options.out.reset(new OutputFile(config.output_file, config.compressor()));
	if (*output_format == Output_format::daa)
		init_daa(*static_cast<OutputFile*>(options.out.get()));
	unique_ptr<OutputFile> unaligned_file, aligned_file;
	if (!config.unaligned.empty())
		unaligned_file = unique_ptr<OutputFile>(new OutputFile(config.unaligned));
	if (!config.aligned_file.empty())
		aligned_file = unique_ptr<OutputFile>(new OutputFile(config.aligned_file));
	timer.finish();

	for (;; ++current_query_chunk) {
		task_timer timer("Loading query sequences", true);

		if (options.self) {
			db_file->set_seqinfo_ptr(query_file_offset);
			options.query.reset(db_file->load_seqs((size_t)(config.chunk_size * 1e9), true, options.db_filter.get()));
			query_file_offset = db_file->tell_seq();
		}
		else
			options.query.reset(new Block(options.query_file->begin(), options.query_file->end(), *format_n, (size_t)(config.chunk_size * 1e9), input_value_traits, config.store_query_quality, false, paired_mode ? 2 : 1));

		if (options.query->empty())
			break;
		timer.finish();
		options.query->seqs().print_stats();
		if ((config.mp_query_chunk >= 0) && (current_query_chunk != (unsigned)config.mp_query_chunk))
			continue;

		if (current_query_chunk == 0 && *output_format != Output_format::daa)
			output_format->print_header(*options.out, align_mode.mode, config.matrix.c_str(), score_matrix.gap_open(), score_matrix.gap_extend(), config.max_evalue, options.query->ids()[0],
				unsigned(align_mode.query_translated ? options.query->source_seqs()[0].length() : options.query->seqs()[0].length()));

		if (config.masking == 1 && !options.self) {
			timer.go("Masking queries");
			mask_seqs(options.query->seqs(), Masking::get());
			timer.finish();
		}

		run_query_chunk(current_query_chunk, *options.out, unaligned_file.get(), aligned_file.get(), options);

		if (file_exists("stop")) {
			message_stream << "Encountered \'stop\' file, shutting down run" << endl;
			break;
		}
	}

	if (options.query_file.unique()) {
		timer.go("Closing the input file");
		for(TextInputFile& f : *options.query_file)
			f.close();
	}

	timer.go("Closing the output file");
	if (*output_format == Output_format::daa) {
		db_file->init_random_access(current_query_chunk, 0);
		finish_daa(*static_cast<OutputFile*>(options.out.get()), *db_file);
		db_file->end_random_access();
	}
	else
		output_format->print_footer(*options.out);
	options.out->finalize();
	if (unaligned_file.get())
		unaligned_file->close();
	if (aligned_file.get())
		aligned_file->close();

	timer.go("Cleaning up");
	options.free();

	timer.finish();
	log_rss();
	message_stream << "Total time = " << total_timer.get() << "s" << endl;
	statistics.print();
	//print_warnings();
}

void run(const shared_ptr<SequenceFile>& db, const shared_ptr<std::list<TextInputFile>>& query, const shared_ptr<Consumer>& out, const shared_ptr<BitVector>& db_filter)
{
	task_timer total;

	align_mode = Align_mode(Align_mode::from_command(config.command));

	message_stream << "Temporary directory: " << TempFile::get_temp_dir() << endl;

	if (config.sensitivity >= Sensitivity::VERY_SENSITIVE)
		::Config::set_option(config.chunk_size, 0.4);
	else
		::Config::set_option(config.chunk_size, 2.0);

	init_output();

	const bool taxon_filter = !config.taxonlist.empty() || !config.taxon_exclude.empty();
	const bool taxon_culling = config.taxon_k != 0;
	SequenceFile::Metadata metadata_flags = SequenceFile::Metadata();
	if (output_format->needs_taxon_id_lists || taxon_filter || taxon_culling)
		metadata_flags |= SequenceFile::Metadata::TAXON_MAPPING;
	if (output_format->needs_taxon_nodes || taxon_filter || taxon_culling)
		metadata_flags |= SequenceFile::Metadata::TAXON_NODES;
	if (output_format->needs_taxon_scientific_names)
		metadata_flags |= SequenceFile::Metadata::TAXON_SCIENTIFIC_NAMES;
	if (output_format->needs_taxon_ranks || taxon_culling)
		metadata_flags |= SequenceFile::Metadata::TAXON_RANKS;

	Config cfg;
	task_timer timer("Opening the database", 1);
	SequenceFile::Flags flags(SequenceFile::Flags::NONE);
	if (flag_any(output_format->flags, Output::Flags::ALL_SEQIDS))
		flags |= SequenceFile::Flags::ALL_SEQIDS;
	if (flag_any(output_format->flags, Output::Flags::FULL_TITLES))
		flags |= SequenceFile::Flags::FULL_TITLES;
	if (flag_any(output_format->flags, Output::Flags::TARGET_SEQS))
		flags |= SequenceFile::Flags::TARGET_SEQS;
	if (db) {
		cfg.db = db;
		if (!query)
			cfg.self = true;
	}
	else
		cfg.db.reset(SequenceFile::auto_create(flags, metadata_flags));
	cfg.query_file = query;
	cfg.db_filter = db_filter;
	cfg.out = out;
	timer.finish();

	message_stream << "Database: " << config.database << ' ';
	message_stream << "(type: " << to_string(cfg.db->type()) << ", ";
	message_stream << "sequences: " << cfg.db->sequence_count() << ", ";
	message_stream << "letters: " << cfg.db->letters() << ')' << endl;
	message_stream << "Block size = " << (size_t)(config.chunk_size * 1e9) << endl;
	score_matrix.set_db_letters(config.db_size ? config.db_size : cfg.db->letters());

	if (output_format->needs_taxon_nodes || taxon_filter || taxon_culling) {
		timer.go("Loading taxonomy nodes");
		cfg.taxon_nodes = cfg.db->taxon_nodes();
		if (taxon_filter) {
			timer.go("Building taxonomy filter");
			cfg.db_filter.reset(cfg.db->filter_by_taxonomy(config.taxonlist, config.taxon_exclude, *cfg.taxon_nodes));
		}
		timer.finish();
	}
	if (output_format->needs_taxon_scientific_names) {
		timer.go("Loading taxonomy names");
		cfg.taxonomy_scientific_names = cfg.db->taxon_scientific_names();
		timer.finish();
	}

	if (!config.seqidlist.empty()) {
		message_stream << "Filtering database by accession list: " << config.seqidlist << endl;
		timer.go("Building database filter");
		cfg.db_filter.reset(cfg.db->filter_by_accession(config.seqidlist));
		timer.finish();
	}

	master_thread(total, cfg);
}

}

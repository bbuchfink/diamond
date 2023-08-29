/****
DIAMOND protein aligner
Copyright (C) 2013-2023 Max Planck Society for the Advancement of Science e.V.
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
#include "../output/output_format.h"
#include "../data/frequent_seeds.h"
#include "../output/daa/daa_write.h"
#include "../data/taxonomy.h"
#include "../masking/masking.h"
#include "../data/block/block.h"
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
#include "../util/async_buffer.h"
#include "config.h"
#include "../data/seed_array.h"
#include "../data/fasta/fasta_file.h"
#ifdef WITH_DNA
#include "../dna/dna_index.h"
#endif


using std::unique_ptr;
using std::endl;
using std::list;
using std::shared_ptr;

namespace Search {

static const int64_t MAX_INDEX_QUERY_SIZE = 32 * MEGABYTES;
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

static string get_ref_part_file_name(const string & prefix, size_t query, string suffix="") {
	if (suffix.size() > 0)
		suffix.append("_");
	const string file_name = append_label(prefix + "_" + suffix, query);
	return join_path(config.parallel_tmpdir, file_name);
}

static string get_ref_block_tmpfile_name(size_t query, size_t block) {
	const string file_name = append_label("ref_block_", query) + append_label("_", block);
	return join_path(config.parallel_tmpdir, file_name);
}

static void run_ref_chunk(SequenceFile &db_file,
	const unsigned query_iteration,
	Consumer &master_out,
	PtrVector<TempFile> &tmp_file,
	Config& cfg)
{
	TaskTimer timer;
	log_rss();
	auto& query_seqs = cfg.query->seqs();

	if (cfg.lin_stage1_target && !config.kmer_ranking && cfg.target.unique()) {
		timer.go("Length sorting reference");
		cfg.target.reset(cfg.target->length_sorted(config.threads_));
	}

	if (config.comp_based_stats == Stats::CBS::COMP_BASED_STATS_AND_MATRIX_ADJUST || flag_any(cfg.output_format->flags, Output::Flags::TARGET_SEQS)) {
		cfg.target->unmasked_seqs() = cfg.target->seqs();
		cfg.target->unmasked_seqs().convert_all_to_std_alph(config.threads_);
	}

	if (cfg.target_masking != MaskingAlgo::NONE && !cfg.lazy_masking) {
		timer.go("Masking reference");
		size_t n = mask_seqs(cfg.target->seqs(), Masking::get(), true, cfg.target_masking);
		timer.finish();
		log_stream << "Masked letters: " << n << endl;
	}

	if (flag_any(cfg.output_format->flags, Output::Flags::SELF_ALN_SCORES)) {
		timer.go("Computing self alignment scores");
		cfg.target->compute_self_aln();
	}

	const bool daa = *cfg.output_format == OutputFormat::daa;
	const bool persist_dict = daa || cfg.iterated();
	if(((cfg.blocked_processing || daa) && !config.global_ranking_targets) || cfg.iterated()) {
		timer.go("Initializing dictionary");
		if (config.multiprocessing || (cfg.current_ref_block == 0 && (!daa || cfg.current_query_block == 0) && query_iteration == 0))
			db_file.init_dict(cfg.current_query_block, cfg.current_ref_block);
		if(!config.global_ranking_targets)
			db_file.init_dict_block(cfg.current_ref_block, cfg.target->seqs().size(), persist_dict);
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
		if (query_seeds_bitset.get()) {
			EnumCfg enum_cfg{ nullptr, 0, 0, cfg.seed_encoding, nullptr, false, false, cfg.seed_complexity_cut, MaskingAlgo::NONE, cfg.minimizer_window, false, false };
			cfg.target->hst() = SeedHistogram(*cfg.target, true, query_seeds_bitset.get(), enum_cfg);
		}
		else if (query_seeds_hashed.get()) {
			EnumCfg enum_cfg{ nullptr, 0, 0, cfg.seed_encoding, nullptr, false, false, cfg.seed_complexity_cut, MaskingAlgo::NONE, cfg.minimizer_window, false, false };
			cfg.target->hst() = SeedHistogram(*cfg.target, true, query_seeds_hashed.get(), enum_cfg);
		}
		else {
			EnumCfg enum_cfg{ nullptr, 0, 0, cfg.seed_encoding, nullptr, false, false, cfg.seed_complexity_cut, cfg.soft_masking, cfg.minimizer_window, false, false };
			cfg.target->hst() = SeedHistogram(*cfg.target, false, &no_filter, enum_cfg);
		}

		timer.go("Allocating buffers");
		char* ref_buffer = SeedArray::alloc_buffer(cfg.target->hst(), cfg.index_chunks),
			* query_buffer = config.target_indexed ? nullptr : SeedArray::alloc_buffer(cfg.query->hst(), cfg.index_chunks);
		timer.finish();

		::HashedSeedSet* target_seeds = nullptr;
		if (config.target_indexed) {
			timer.go("Loading database seed index");
			target_seeds = new ::HashedSeedSet(db_file.file_name() + ".seed_idx");
			timer.finish();
		}
        if((config.command != ::Config::blastn)){
            for (unsigned i = 0; i < shapes.count(); ++i) {
                if(config.global_ranking_targets)
                    cfg.global_ranking_buffer.reset(new Config::RankingBuffer());
                search_shape(i, cfg.current_query_block, query_iteration, query_buffer, ref_buffer, cfg, target_seeds); //index_targets(0,cfg,ref_buffer,target_seeds);
                if (config.global_ranking_targets)
                    Extension::GlobalRanking::update_table(cfg);
            }
        }
#ifdef WITH_DNA
        else
            cfg.dna_ref_index = std::make_unique<Dna::Index>(cfg, ref_buffer);
#endif

        log_rss();
        timer.go("Deallocating buffers");
#ifdef WITH_DNA
        if(config.command != ::Config::blastn)
#endif
		delete[] ref_buffer;
		delete[] query_buffer;
		delete target_seeds;

		timer.go("Clearing query masking");
		Frequent_seeds::clear_masking(query_seqs);
		timer.finish();
		log_rss();
	}

    Consumer* out;
	const bool temp_output = (cfg.blocked_processing || cfg.iterated()) && !config.global_ranking_targets;
	if (temp_output) {
		timer.go("Opening temporary output file");
		if (config.multiprocessing) {
			const string file_name = get_ref_block_tmpfile_name(cfg.current_query_block, cfg.current_ref_block);
			tmp_file.push_back(new TempFile(file_name));
		} else {
			tmp_file.push_back(new TempFile());
		}
		out = &tmp_file.back();
	}
	else
		out = &master_out;

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

static void run_query_iteration(const unsigned query_iteration,
	Consumer& master_out,
	OutputFile* unaligned_file,
	OutputFile* aligned_file,
	PtrVector<TempFile>& tmp_file,
	Config& options)
{
	TaskTimer timer;
	auto P = Parallelizer::get();
	auto& query_seqs = options.query->seqs();
	auto& db_file = *options.db;
	if (query_iteration > 0)
		options.query_skip.reset(new vector<bool> (query_aligned));

	if (config.algo == ::Config::Algo::AUTO &&
		(!sensitivity_traits[(int)align_mode.sequence_type].at(config.sensitivity).support_query_indexed
			|| query_seqs.letters() > MAX_INDEX_QUERY_SIZE
			|| options.db_letters < MIN_QUERY_INDEXED_DB_SIZE
			|| config.target_indexed
			|| config.swipe_all
			|| options.minimizer_window))
		config.algo = ::Config::Algo::DOUBLE_INDEXED;
	if (config.algo == ::Config::Algo::AUTO || config.algo == ::Config::Algo::QUERY_INDEXED) {
		timer.go("Building query seed set");
		query_seeds_hashed.reset(new HashedSeedSet(*options.query, options.query_skip.get(), options.seed_complexity_cut, options.soft_masking));
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
		query_seeds_bitset.reset(new SeedSet(*options.query, 1.0, options.query_skip.get(), options.seed_complexity_cut, options.soft_masking));
		options.seed_encoding = SeedEncoding::CONTIGUOUS;
		timer.finish();
	}

	const Sensitivity sens = options.sensitivity[query_iteration].sensitivity;
	::Config::set_option(options.index_chunks, config.lowmem_, 0u, config.algo == ::Config::Algo::DOUBLE_INDEXED ? sensitivity_traits[(int)align_mode.sequence_type].at(sens).index_chunks : 1u);
	options.lazy_masking = config.algo != ::Config::Algo::DOUBLE_INDEXED && options.target_masking != MaskingAlgo::NONE && config.frame_shift == 0;
	if (config.command != ::Config::blastn && options.gapped_filter_evalue != 0.0) {
		options.cutoff_gapped1 = { config.gapped_filter_evalue1 };
		options.cutoff_gapped2 = { options.gapped_filter_evalue };
		options.cutoff_gapped1_new = { config.gapped_filter_evalue1 };
		options.cutoff_gapped2_new = { options.gapped_filter_evalue };
	}

	if (options.current_query_block == 0 && query_iteration == 0) {
		message_stream << "Algorithm: " << to_string(config.algo) << endl;
		if (config.freq_masking && !config.lin_stage1)
			verbose_stream << "Seed frequency SD: " << options.freq_sd << endl;
		verbose_stream << "Shape configuration: " << ::shapes << endl;
	}	

	if (config.global_ranking_targets) {
		timer.go("Allocating global ranking table");
		options.ranking_table.reset(new Search::Config::RankingTable(query_seqs.size() * config.global_ranking_targets / align_mode.query_contexts));
	}

	if (!config.swipe_all && !config.target_indexed) {
		timer.go("Building query histograms");
		EnumCfg enum_cfg{ nullptr, 0, 0, options.seed_encoding, options.query_skip.get(), false, false, options.seed_complexity_cut,
			options.soft_masking, options.minimizer_window, (bool)query_seeds_hashed, false };
		options.query->hst() = SeedHistogram(*options.query, false, &no_filter, enum_cfg);
		timer.finish();
	}

	log_rss();
	bool mp_last_chunk = false;
	//const bool lazy_masking = config.algo == ::Config::Algo::QUERY_INDEXED && (config.masking == 1 || config.target_seg == 1) && (config.global_ranking_targets == 0);
	SequenceFile::LoadFlags load_flags = SequenceFile::LoadFlags::SEQS;
	if (!flag_any(db_file.format_flags(), SequenceFile::FormatFlags::TITLES_LAZY) || config.no_self_hits)
		load_flags |= SequenceFile::LoadFlags::TITLES;
	if (options.lazy_masking)
		load_flags |= SequenceFile::LoadFlags::LAZY_MASKING;

	if (config.multiprocessing) {
		db_file.set_seqinfo_ptr(0);
		P->create_stack_from_file(stack_align_todo, get_ref_part_file_name(stack_align_todo, options.current_query_block));
		auto work = P->get_stack(stack_align_todo);
		P->create_stack_from_file(stack_align_wip, get_ref_part_file_name(stack_align_wip, options.current_query_block));
		auto wip = P->get_stack(stack_align_wip);
		P->create_stack_from_file(stack_align_done, get_ref_part_file_name(stack_align_done, options.current_query_block));
		auto done = P->get_stack(stack_align_done);
		P->create_stack_from_file(stack_join_todo, get_ref_part_file_name(stack_join_todo, options.current_query_block));
		auto join_work = P->get_stack(stack_join_todo);

		string buf;

		while ((!file_exists("stop")) && (work->pop(buf))) {
			wip->push(buf);

			Chunk chunk = to_chunk(buf);

			P->log("SEARCH BEGIN " + std::to_string(options.current_query_block) + " " + std::to_string(chunk.i));

			options.target.reset(db_file.load_seqs((size_t)(0), options.db_filter.get(), load_flags, chunk));
			options.current_ref_block = chunk.i;
			options.blocked_processing = true;
			if (!config.mp_self || chunk.i >= options.current_query_block)
				run_ref_chunk(db_file, query_iteration, master_out, tmp_file, options);
			else {
				const string file_name = get_ref_block_tmpfile_name(options.current_query_block, options.current_ref_block);
				tmp_file.push_back(new TempFile(file_name));
				tmp_file.back().write(IntermediateRecord::FINISHED);
				db_file.init_dict(options.current_query_block, options.current_ref_block);
				db_file.close_dict_block(false);
			}

			tmp_file.back().close();

			size_t size_after_push = 0;
			done->push(buf, size_after_push);
			if (size_after_push == (size_t)db_file.get_n_partition_chunks()) {
				join_work->push("TOKEN");
			}
			wip->remove(buf);

			P->log("SEARCH END " + std::to_string(options.current_query_block) + " " + std::to_string(chunk.i));
			log_rss();
		}

		tmp_file.clear();
		P->delete_stack(stack_align_todo);
		P->delete_stack(stack_align_wip);
		P->delete_stack(stack_align_done);
	}
	else {
		/*if (config.self && !config.lin_stage1 && !db_file.eof())
			db_file.set_seqinfo_ptr(options.query->oid_end());
		else if (!config.self || options.current_query_block != 0 || !db_file.eof())
			db_file.set_seqinfo_ptr(0);*/
		timer.go("Seeking in database");
		db_file.set_seqinfo_ptr((config.self && !config.lin_stage1) ? options.query->oid_end() : 0);
		timer.finish();
		for (options.current_ref_block = 0; ; ++options.current_ref_block) {
			if (config.self && ((config.lin_stage1 && options.current_ref_block == options.current_query_block) || (!config.lin_stage1 && options.current_ref_block == 0))) {
				options.target = options.query;
				if (config.lin_stage1) {
					timer.go("Seeking in database");
					db_file.set_seqinfo_ptr(options.query->oid_end());
					timer.finish();
				}
			}
			else {
				timer.go("Loading reference sequences");
				options.target.reset(db_file.load_seqs(config.block_size(), options.db_filter.get(), load_flags));
			}
			if (options.current_ref_block == 0) {
				db_file.reopen();
				const int64_t db_seq_count = options.db_filter ? options.db_filter->one_count() : options.db->sequence_count();
				options.blocked_processing = config.global_ranking_targets || options.target->seqs().size() < db_seq_count;
			}
			if (options.target->empty()) break;
			timer.finish();
			run_ref_chunk(db_file, query_iteration, master_out, tmp_file, options);
		}
		log_rss();
	}

	timer.go("Deallocating buffers");
	query_seeds_hashed.reset();
	query_seeds_bitset.reset();
	options.query_skip.reset();

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

static void run_query_chunk(Consumer &master_out,
	OutputFile *unaligned_file,
	OutputFile *aligned_file,
	Config &options)
{
	auto P = Parallelizer::get();
	TaskTimer timer;
	auto& db_file = *options.db;
	auto& query_seqs = options.query->seqs();

	PtrVector<TempFile> tmp_file;
	if (options.track_aligned_queries) {
		query_aligned.clear();
		query_aligned.insert(query_aligned.end(), options.query->source_seq_count(), false);
	}
	if (flag_any(options.output_format->flags, Output::Flags::SELF_ALN_SCORES)) {
		timer.go("Computing self alignment scores");
		options.query->compute_self_aln();
	}

	log_rss();

	BlockId aligned = 0;
	for (unsigned query_iteration = 0; query_iteration < options.sensitivity.size() && aligned < options.query->source_seq_count(); ++query_iteration) {
		setup_search(options.sensitivity[query_iteration].sensitivity, options);
		options.lin_stage1_target = config.linsearch || options.sensitivity[query_iteration].linearize;
		run_query_iteration(query_iteration, master_out, unaligned_file, aligned_file, tmp_file, options);
		if (options.iterated()) {
			aligned += options.iteration_query_aligned;
			message_stream << "Aligned " << options.iteration_query_aligned << '/' << options.query->source_seq_count() << " queries in this iteration, "
				<< aligned << '/' << options.query->source_seq_count() << " total." << endl;
			options.iteration_query_aligned = 0;
		}
	}

	log_rss();

	if (options.blocked_processing || config.multiprocessing || options.iterated()) {
		if(!config.global_ranking_targets) timer.go("Joining output blocks");

		if (config.multiprocessing) {
			P->create_stack_from_file(stack_join_todo, get_ref_part_file_name(stack_join_todo, options.current_query_block));
			auto work = P->get_stack(stack_join_todo);

			string buf;

			if ((! file_exists("stop")) && (work->pop(buf) > 0)) {
				P->log("JOIN BEGIN "+std::to_string(options.current_query_block));

				P->create_stack_from_file(stack_join_wip, get_ref_part_file_name(stack_join_wip, options.current_query_block));
				auto wip = P->get_stack(stack_join_wip);
				wip->clear();
				P->create_stack_from_file(stack_join_done, get_ref_part_file_name(stack_join_done, options.current_query_block));
				auto done = P->get_stack(stack_join_done);
				done->clear();

				wip->push(buf);
				work->clear();

				options.current_ref_block = db_file.get_n_partition_chunks();

				vector<string> tmp_file_names;
				for (int64_t i=0; i<options.current_ref_block; ++i) {
					tmp_file_names.push_back(get_ref_block_tmpfile_name(options.current_query_block, i));
				}

				const string query_chunk_output_file = append_label(config.output_file + "_", options.current_query_block);
				Consumer *query_chunk_out(new OutputFile(query_chunk_output_file, config.compressor()));
				// if (*output_format != Output_format::daa)
				// 	output_format->print_header(*query_chunk_out, align_mode.mode, config.matrix.c_str(), score_matrix.gap_open(), score_matrix.gap_extend(), config.max_evalue, query_ids::get()[0],
				// 		unsigned(align_mode.query_translated ? query_source_seqs::get()[0].length() : query_seqs::get()[0].length()));

				join_blocks(options.current_ref_block, *query_chunk_out, tmp_file, options, db_file, tmp_file_names);

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

				P->log("JOIN END "+std::to_string(options.current_query_block));
			}
			P->delete_stack(stack_join_todo);
		} else {
			if (!tmp_file.empty())
				join_blocks(options.current_ref_block, master_out, tmp_file, options, db_file);
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

static void master_thread(TaskTimer &total_timer, Config &options)
{
	log_rss();
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
		db_file->create_partition_balanced((size_t)(config.chunk_size*1e9));
	}

	TaskTimer timer("Opening the input file", true);
	if (!options.self) {
		if (config.query_file.empty() && !options.query_file) {
			std::cerr << "Query file parameter (--query/-q) is missing. Input will be read from stdin." << endl;
			config.query_file.push_back("");
		}
		if (!options.query_file)
			options.query_file.reset(new FastaFile(config.query_file, SequenceFile::Metadata(), SequenceFile::Flags(), input_value_traits));
	}

	options.current_query_block = 0;
	OId query_file_offset = 0;
	SequenceFile::LoadFlags load_flags = SequenceFile::LoadFlags::ALL;
	if (config.store_query_quality)
		load_flags |= SequenceFile::LoadFlags::QUALITY;
	if (config.command == ::Config::blastn)
		load_flags |= SequenceFile::LoadFlags::DNA_PRESERVATION;
	if (config.multiprocessing && config.mp_init) {
		TaskTimer timer("Counting query blocks", true);

		size_t block_count = 0;
		do {
			if (options.self) {
				db_file->set_seqinfo_ptr(query_file_offset);
				options.query.reset(db_file->load_seqs((size_t)(config.chunk_size * 1e9), options.db_filter.get(), SequenceFile::LoadFlags::ALL));
				query_file_offset = db_file->tell_seq();
			} else
				options.query.reset(options.query_file->load_seqs((int64_t)(config.chunk_size * 1e9), nullptr, load_flags));
			++block_count;
		} while (!options.query->empty());
		if (options.self) {
			db_file->set_seqinfo_ptr(0);
			query_file_offset = 0;
		} else {
			options.query_file->set_seqinfo_ptr(0);
		}

		for (size_t i = 0; i < block_count - 1; ++i) {
			const string annotation = "# query_chunk=" + std::to_string(i);
			db_file->save_partition(get_ref_part_file_name(stack_align_todo, i), annotation);
		}

		return;
	}

	timer.go("Opening the output file");
	if (!options.out)
		options.out.reset(new OutputFile(config.output_file, config.compressor()));
	if (*options.output_format == OutputFormat::daa)
		init_daa(*static_cast<OutputFile*>(options.out.get()));
	unique_ptr<OutputFile> unaligned_file, aligned_file;
	if (!config.unaligned.empty())
		unaligned_file = unique_ptr<OutputFile>(new OutputFile(config.unaligned));
	if (!config.aligned_file.empty())
		aligned_file = unique_ptr<OutputFile>(new OutputFile(config.aligned_file));
	timer.finish();

	for (;query_file_offset < db_file->sequence_count(); ++options.current_query_block) {
		log_rss();

		if (options.self) {
			timer.go("Seeking in database");
			db_file->set_seqinfo_ptr(query_file_offset);
			timer.finish();
			timer.go("Loading query sequences");
			options.query.reset(db_file->load_seqs(config.block_size(), options.db_filter.get(), SequenceFile::LoadFlags::ALL));
			query_file_offset = db_file->tell_seq();
		}
		else {
			timer.go("Loading query sequences");
			options.query.reset(options.query_file->load_seqs(config.block_size(), nullptr, load_flags));
		}
		timer.finish();

		if (options.query->empty())
			break;
		options.query->seqs().print_stats();
		if ((config.mp_query_chunk >= 0) && (options.current_query_block != config.mp_query_chunk))
			continue;

#ifndef KEEP_TARGET_ID
		if (config.lin_stage1 && !config.kmer_ranking) {
			timer.go("Length sorting queries");
			options.query.reset(options.query->length_sorted(config.threads_));
			timer.finish();
		}
#endif

		if (options.current_query_block == 0 && *options.output_format != OutputFormat::daa && options.query->has_ids())
			options.output_format->print_header(*options.out, align_mode.mode, config.matrix.c_str(), score_matrix.gap_open(), score_matrix.gap_extend(), config.max_evalue, options.query->ids()[0],
				unsigned(align_mode.query_translated ? options.query->source_seqs()[0].length() : options.query->seqs()[0].length()));

		if (options.query_masking != MaskingAlgo::NONE) {
			timer.go("Masking queries");
			mask_seqs(options.query->seqs(), Masking::get(), true, options.query_masking);
			timer.finish();
		}

		run_query_chunk(*options.out, unaligned_file.get(), aligned_file.get(), options);

		if (file_exists("stop")) {
			message_stream << "Encountered \'stop\' file, shutting down run" << endl;
			break;
		}
	}

	if (options.query_file.unique()) {
		timer.go("Closing the input file");
		options.query_file->close();
	}

	timer.go("Closing the output file");
	if (*options.output_format == OutputFormat::daa) {
		db_file->init_random_access(options.current_query_block, 0);
		finish_daa(*static_cast<OutputFile*>(options.out.get()), *db_file);
		db_file->end_random_access();
	}
	else
		options.output_format->print_footer(*options.out);
	options.out->finalize();
	if (unaligned_file.get())
		unaligned_file->close();
	if (aligned_file.get())
		aligned_file->close();

	if (!config.unaligned_targets.empty()) {
		timer.go("Writing unaligned targets");
		options.db->write_accession_list(options.aligned_targets, config.unaligned_targets);
	}

	timer.go("Closing the database");
	options.db.reset();

	timer.go("Cleaning up");
	options.free();

	timer.finish();
	log_rss();
	message_stream << "Total time = " << total_timer.get() << "s" << endl;
	statistics.print();
	//print_warnings();
}

void run(const shared_ptr<SequenceFile>& db, const shared_ptr<SequenceFile>& query, const shared_ptr<Consumer>& out, const shared_ptr<BitVector>& db_filter)
{
	TaskTimer total;

	align_mode = AlignMode(AlignMode::from_command(config.command));
    (align_mode.sequence_type == SequenceType::amino_acid) ? value_traits = amino_acid_traits : value_traits = nucleotide_traits;

	message_stream << "Temporary directory: " << TempFile::get_temp_dir() << endl;

	if (config.sensitivity >= Sensitivity::VERY_SENSITIVE)
		::Config::set_option(config.chunk_size, 0.4);
	else
		::Config::set_option(config.chunk_size, 2.0);

	Config cfg;
	cfg.output_format.reset(init_output(cfg.max_target_seqs));
	statistics.reset();

	const bool taxon_filter = !config.taxonlist.empty() || !config.taxon_exclude.empty();
	const bool taxon_culling = config.taxon_k != 0;
	SequenceFile::Metadata metadata_flags = SequenceFile::Metadata();
	if (cfg.output_format->needs_taxon_id_lists || taxon_filter || taxon_culling)
		metadata_flags |= SequenceFile::Metadata::TAXON_MAPPING;
	if (cfg.output_format->needs_taxon_nodes || taxon_filter || taxon_culling)
		metadata_flags |= SequenceFile::Metadata::TAXON_NODES;
	if (cfg.output_format->needs_taxon_scientific_names)
		metadata_flags |= SequenceFile::Metadata::TAXON_SCIENTIFIC_NAMES;
	if (cfg.output_format->needs_taxon_ranks || taxon_culling)
		metadata_flags |= SequenceFile::Metadata::TAXON_RANKS;

	TaskTimer timer;
	SequenceFile::Flags flags(SequenceFile::Flags::NEED_LETTER_COUNT);
	if (flag_any(cfg.output_format->flags, Output::Flags::ALL_SEQIDS))
		flags |= SequenceFile::Flags::ALL_SEQIDS;
	if (flag_any(cfg.output_format->flags, Output::Flags::FULL_TITLES) || config.no_self_hits)
		flags |= SequenceFile::Flags::FULL_TITLES;
	if (flag_any(cfg.output_format->flags, Output::Flags::TARGET_SEQS))
		flags |= SequenceFile::Flags::TARGET_SEQS;
	if (flag_any(cfg.output_format->flags, Output::Flags::SELF_ALN_SCORES))
		flags |= SequenceFile::Flags::SELF_ALN_SCORES;
	if (!config.unaligned_targets.empty())
		flags |= SequenceFile::Flags::OID_TO_ACC_MAPPING;
	if (db) {
		cfg.db = db;
		if (!query)
			cfg.self = true;
	}
	else {
		timer.go("Opening the database");
		cfg.db.reset(SequenceFile::auto_create({ config.database }, flags, metadata_flags, value_traits));
	}
	if (config.multiprocessing && cfg.db->type() == SequenceFile::Type::FASTA)
		throw std::runtime_error("Multiprocessing mode is not compatible with FASTA databases.");
	cfg.db_seqs = cfg.db->sequence_count();
	cfg.db_letters = cfg.db->letters();
	cfg.ref_blocks = cfg.db->total_blocks();
	cfg.query_file = query;
	cfg.db_filter = db_filter;
	cfg.out = out;
	if (!config.unaligned_targets.empty())
		cfg.aligned_targets.insert(cfg.aligned_targets.begin(), cfg.db->sequence_count(), false);
	timer.finish();

	message_stream << "Database: " << config.database << ' ';
	message_stream << "(type: " << to_string(cfg.db->type()) << ", ";
	message_stream << "sequences: " << cfg.db->sequence_count() << ", ";
	message_stream << "letters: " << cfg.db->letters() << ')' << endl;
	message_stream << "Block size = " << (size_t)(config.chunk_size * 1e9) << endl;
	score_matrix.set_db_letters(config.db_size ? config.db_size : cfg.db->letters());

	if (cfg.output_format->needs_taxon_nodes || taxon_filter || taxon_culling) {
		if (taxon_filter) {
			timer.go("Building taxonomy filter");
			cfg.db_filter.reset(cfg.db->filter_by_taxonomy(config.taxonlist, config.taxon_exclude));
		}
		timer.finish();
	}

	if (!config.seqidlist.empty()) {
		if (taxon_filter)
			throw std::runtime_error("--seqidlist is not compatible with taxonomy filtering.");
		message_stream << "Filtering database by accession list: " << config.seqidlist << endl;
		timer.go("Building database filter");
		cfg.db_filter.reset(cfg.db->filter_by_accession(config.seqidlist));
		timer.finish();
	}

	if (align_mode.sequence_type == SequenceType::nucleotide)
		cfg.score_builder.reset(new Stats::Blastn_Score(config.match_reward, config.mismatch_penalty, config.gap_open, config.gap_extend, cfg.db_letters, cfg.db->sequence_count()));


    master_thread(total, cfg);
	log_rss();
}

}

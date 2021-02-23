/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
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
#include "../output/daa_write.h"
#include "../data/taxonomy.h"
#include "../basic/masking.h"
#include "../data/ref_dictionary.h"
#include "../data/metadata.h"
#include "../search/search.h"
#include "workflow.h"
#include "../util/io/consumer.h"
#include "../util/parallel/thread_pool.h"
#include "../util/parallel/multiprocessing.h"
#include "../util/parallel/parallelizer.h"
#include "../util/system/system.h"
#include "../align/target.h"
#include "../data/enum_seeds.h"

using std::unique_ptr;
using std::endl;

namespace Workflow { namespace Search {

static const string label_align = "align";
static const string stack_align_todo = label_align + "_todo";
static const string stack_align_wip = label_align + "_wip";
static const string stack_align_done = label_align + "_done";

static const string label_join = "join";
static const string stack_join_todo = label_join + "_todo";
static const string stack_join_wip = label_join + "_wip";
static const string stack_join_done = label_join + "_done";


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
	unsigned query_chunk,
	pair<size_t, size_t> query_len_bounds,
	char *query_buffer,
	Consumer &master_out,
	PtrVector<TempFile> &tmp_file,
	const Parameters &params,
	const Metadata &metadata)
{
	log_rss();

	if (config.comp_based_stats == Stats::CBS::COMP_BASED_STATS_AND_MATRIX_ADJUST)
		ref_seqs_unmasked::data_ = new SequenceSet(*ref_seqs::data_);

	task_timer timer;
	if (config.masking == 1 && !config.no_ref_masking) {
		timer.go("Masking reference");
		size_t n = mask_seqs(*ref_seqs::data_, Masking::get());
		timer.finish();
		log_stream << "Masked letters: " << n << endl;
	}

	ReferenceDictionary::get().init(safe_cast<unsigned>(ref_seqs::get().get_length()), block_to_database_id);

	timer.go("Initializing temporary storage");
	Trace_pt_buffer::instance = new Trace_pt_buffer(query_seqs::data_->get_length() / align_mode.query_contexts,
		config.tmpdir,
		config.query_bins);

	if (!config.swipe_all) {
		timer.go("Building reference histograms");
		if (config.algo == Config::query_indexed)
			ref_hst = Partitioned_histogram(*ref_seqs::data_, false, query_seeds);
		else if (query_seeds_hashed != 0)
			ref_hst = Partitioned_histogram(*ref_seqs::data_, true, query_seeds_hashed);
		else
			ref_hst = Partitioned_histogram(*ref_seqs::data_, false, &no_filter);

		timer.go("Allocating buffers");
		char *ref_buffer = SeedArray::alloc_buffer(ref_hst);
		timer.finish();

		HashedSeedSet* target_seeds = nullptr;
		if (config.target_indexed) {
			timer.go("Loading database seed index");
			target_seeds = new HashedSeedSet(config.database + ".seed_idx");
			timer.finish();
		}

		for (unsigned i = 0; i < shapes.count(); ++i)
			search_shape(i, query_chunk, query_buffer, ref_buffer, params, target_seeds);

		timer.go("Deallocating buffers");
		delete[] ref_buffer;
		delete target_seeds;

		timer.go("Clearing query masking");
		Frequent_seeds::clear_masking(*query_seqs::data_);
	}

	Consumer* out;
	if (blocked_processing) {
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

	if (config.target_seg == 1) {
		timer.go("SEG masking targets");
		mask_seqs(*ref_seqs::data_, Masking::get(), true, Masking::Algo::SEG);
	}

	timer.go("Computing alignments");
	align_queries(*Trace_pt_buffer::instance, out, params, metadata);
	delete Trace_pt_buffer::instance;

	if (blocked_processing)
		IntermediateRecord::finish_file(*out);

	timer.go("Deallocating reference");
	delete ref_seqs::data_;
	delete ref_ids::data_;
	if (config.comp_based_stats == Stats::CBS::COMP_BASED_STATS_AND_MATRIX_ADJUST)
		delete ref_seqs_unmasked::data_;
	timer.finish();
}

void deallocate_queries() {
	delete query_seqs::data_;
	delete query_ids::data_;
	delete query_source_seqs::data_;
	delete query_qual;
}

void run_query_chunk(SequenceFile &db_file,
	unsigned query_chunk,
	Consumer &master_out,
	OutputFile *unaligned_file,
	OutputFile *aligned_file,
	const Metadata &metadata,
	const Options &options)
{
	auto P = Parallelizer::get();

	task_timer timer("Building query seed set");
	if (query_chunk == 0)
		setup_search_cont();
	if (config.algo == -1) {
		if (config.sensitivity >= Sensitivity::VERY_SENSITIVE || config.sensitivity == Sensitivity::MID_SENSITIVE || config.sensitivity == Sensitivity::FAST) {
			config.algo = Config::double_indexed;
		}
		else {
			query_seeds = new Seed_set(query_seqs::get(), SINGLE_INDEXED_SEED_SPACE_MAX_COVERAGE);
			timer.finish();
			log_stream << "Seed space coverage = " << query_seeds->coverage() << endl;
			if (use_single_indexed(query_seeds->coverage(), query_seqs::get().letters(), db_file.letters()))
				config.algo = Config::query_indexed;
			else {
				config.algo = Config::double_indexed;
				delete query_seeds;
				query_seeds = NULL;
			}
		}
	}
	else if (config.algo == Config::query_indexed) {
		if (config.sensitivity >= Sensitivity::VERY_SENSITIVE)
			throw std::runtime_error("Query-indexed algorithm not available for this sensitivity setting.");
		query_seeds = new Seed_set(query_seqs::get(), 2);
		timer.finish();
		log_stream << "Seed space coverage = " << query_seeds->coverage() << endl;
	}
	else
		timer.finish();
	if (query_chunk == 0)
		setup_search();
	if (config.algo == Config::double_indexed && config.small_query) {
		timer.go("Building query seed hash set");
		query_seeds_hashed = new HashedSeedSet(query_seqs::get());
	}
	timer.finish();

	const Parameters params{
	db_file.sequence_count(),
	db_file.letters(),
	db_file.total_blocks(),
	config.gapped_filter_evalue1,
	config.gapped_filter_evalue,
	config.gapped_filter_evalue1,
	config.gapped_filter_evalue
	};

	char *query_buffer = nullptr;
	const pair<size_t, size_t> query_len_bounds = query_seqs::data_->len_bounds(shapes[0].length_);

	if (!config.swipe_all && !config.target_indexed) {
		timer.go("Building query histograms");
		query_hst = Partitioned_histogram(*query_seqs::data_, false, &no_filter);

		timer.go("Allocating buffers");
		query_buffer = SeedArray::alloc_buffer(query_hst);
		timer.finish();
	}

	log_rss();

	PtrVector<TempFile> tmp_file;
	query_aligned.clear();
	query_aligned.insert(query_aligned.end(), query_ids::get().get_length(), false);
	if(config.query_memory)
		Extension::memory = new Extension::Memory(query_ids::get().get_length());
	db_file.set_seqinfo_ptr(0);
	Chunk chunk;
	bool mp_last_chunk = false;

	log_rss();

	if (config.multiprocessing) {
		auto work = P->get_stack(stack_align_todo);
		P->create_stack_from_file(stack_align_wip, get_ref_part_file_name(stack_align_wip, query_chunk));
		auto wip = P->get_stack(stack_align_wip);
		P->create_stack_from_file(stack_align_done, get_ref_part_file_name(stack_align_done, query_chunk));
		auto done = P->get_stack(stack_align_done);
		string buf;

		while (work->pop(buf)) {
			wip->push(buf);

			Chunk chunk = to_chunk(buf);

			P->log("SEARCH BEGIN "+std::to_string(query_chunk)+" "+std::to_string(chunk.i));

			db_file.load_seqs(&block_to_database_id, (size_t)(0), &ref_seqs::data_, &ref_ids::data_, true, options.db_filter ? options.db_filter : metadata.taxon_filter, true, chunk);
			run_ref_chunk(db_file, query_chunk, query_len_bounds, query_buffer, master_out, tmp_file, params, metadata);

			ReferenceDictionary::get().save_block(query_chunk, chunk.i);
			ReferenceDictionary::get().clear_block(chunk.i);

			size_t size_after_push = 0;
			done->push(buf, size_after_push);
			if (size_after_push == db_file.get_n_partition_chunks()) {
				mp_last_chunk = true;
			}
			wip->pop(buf);

			P->log("SEARCH END "+std::to_string(query_chunk)+" "+std::to_string(chunk.i));
			log_rss();
		}
	} else {
		for (current_ref_block = 0;
			 db_file.load_seqs(&block_to_database_id, (size_t)(config.chunk_size*1e9), &ref_seqs::data_, &ref_ids::data_, true, options.db_filter ? options.db_filter : metadata.taxon_filter);
			 ++current_ref_block) {
			run_ref_chunk(db_file, query_chunk, query_len_bounds, query_buffer, master_out, tmp_file, params, metadata);
		}
		log_rss();
	}

	timer.go("Deallocating buffers");
	delete[] query_buffer;
	delete query_seeds;
	delete Extension::memory;
	query_seeds = 0;

	log_rss();

	if (config.multiprocessing) {
		for (auto f : tmp_file) {
			f->close();
		}
		tmp_file.clear();
	}

	log_rss();

	if (blocked_processing) {
		timer.go("Joining output blocks");

		if (config.multiprocessing) {

			if (mp_last_chunk) {
				P->create_stack_from_file(stack_join_todo, get_ref_part_file_name(stack_join_todo, query_chunk));
				auto work = P->get_stack(stack_join_todo);

				P->create_stack_from_file(stack_join_wip, get_ref_part_file_name(stack_join_wip, query_chunk));
				auto wip = P->get_stack(stack_join_wip);

				P->create_stack_from_file(stack_join_done, get_ref_part_file_name(stack_join_done, query_chunk));
				auto done = P->get_stack(stack_join_done);
				string buf;

				work->pop(buf);
				wip->push(buf);
				P->log("JOIN BEGIN "+std::to_string(query_chunk));

				current_ref_block = db_file.get_n_partition_chunks();

				ReferenceDictionary::get().restore_blocks(query_chunk, current_ref_block);

				vector<string> tmp_file_names;
				for (size_t i=0; i<current_ref_block; ++i) {
					tmp_file_names.push_back(get_ref_block_tmpfile_name(query_chunk, i));
				}

				const string query_chunk_output_file = append_label(config.output_file + "_", current_query_chunk);
				Consumer *query_chunk_out(new OutputFile(query_chunk_output_file, config.compression == 1));
				// if (*output_format != Output_format::daa)
				// 	output_format->print_header(*query_chunk_out, align_mode.mode, config.matrix.c_str(), score_matrix.gap_open(), score_matrix.gap_extend(), config.max_evalue, query_ids::get()[0],
				// 		unsigned(align_mode.query_translated ? query_source_seqs::get()[0].length() : query_seqs::get()[0].length()));

				join_blocks(current_ref_block, *query_chunk_out, tmp_file, params, metadata, db_file, tmp_file_names);

				// if (*output_format == Output_format::daa)
				// 	// finish_daa(*static_cast<OutputFile*>(query_chunk_out), *db_file);
				// 	throw std::runtime_error("output_format::daa");
				// else
				// 	output_format->print_footer(*query_chunk_out);
				query_chunk_out->finalize();
				delete query_chunk_out;

				ReferenceDictionary::get().clear_block_instances();
				for (auto f : tmp_file_names) {
					std::remove(f.c_str());
				}

				done->push(buf);
				wip->pop(buf);
				P->log("JOIN END "+std::to_string(query_chunk));
			}
			P->delete_stack(stack_align_done);

		} else {
			join_blocks(current_ref_block, master_out, tmp_file, params, metadata, db_file);
		}
	}

	if (unaligned_file) {
		timer.go("Writing unaligned queries");
		write_unaligned(unaligned_file);
	}
	if (aligned_file) {
		timer.go("Writing aligned queries");
		write_aligned(aligned_file);
	}

	timer.go("Deallocating queries");
	deallocate_queries();
	if (*output_format != Output_format::daa)
		ReferenceDictionary::get().clear();
}



void master_thread(SequenceFile *db_file, task_timer &total_timer, Metadata &metadata, const Options &options)
{
	auto P = Parallelizer::get();
	if (config.multiprocessing) {
		P->init(config.parallel_tmpdir);
		db_file->create_partition_balanced((size_t)(config.chunk_size*1e9));
	}

	task_timer timer("Opening the input file", true);
	list<TextInputFile> *query_file = nullptr;
	const Sequence_file_format *format_n = nullptr;
	bool paired_mode = false;
	if (!options.self) {
		if (config.query_file.empty() && !options.query_file)
			std::cerr << "Query file parameter (--query/-q) is missing. Input will be read from stdin." << endl;
		if (options.query_file)
			query_file = options.query_file;
		else {
			query_file = new list<TextInputFile>;
			for(const string& f : config.query_file)
				query_file->emplace_back(f);
			if (query_file->empty())
				query_file->emplace_back("");
			paired_mode = query_file->size() == 2;
		}
		format_n = guess_format(query_file->front());
	}

	current_query_chunk = 0;
	size_t query_file_offset = 0;

	if (config.multiprocessing && config.mp_init) {
		task_timer timer("Counting query blocks", true);

		for (;; ++current_query_chunk) {
			if (options.self) {
				db_file->set_seqinfo_ptr(query_file_offset);
				if (!db_file->load_seqs(&query_block_to_database_id, (size_t)(config.chunk_size * 1e9),
					&query_seqs::data_,
					&query_ids::data_,
					true,
					options.db_filter))
					break;
				deallocate_queries();
				query_file_offset = db_file->tell_seq();
			} else {
				if (!load_seqs(query_file->begin(), query_file->end(), *format_n, &query_seqs::data_, query_ids::data_, &query_source_seqs::data_,
					config.store_query_quality ? &query_qual : nullptr,
					(size_t)(config.chunk_size * 1e9), config.qfilt, input_value_traits, paired_mode ? 2 : 1))
					break;
				deallocate_queries();
			}
		}
		if (options.self) {
			db_file->set_seqinfo_ptr(0);
			query_file_offset = 0;
		} else {
			query_file->front().rewind();
		}

		for (size_t i = 0; i < current_query_chunk; ++i) {
			const string annotation = "# query_chunk=" + std::to_string(i);
			db_file->save_partition(get_ref_part_file_name(stack_align_todo, i), annotation);

			P->create_stack_from_file(stack_join_todo, get_ref_part_file_name(stack_join_todo, i));
			P->get_stack(stack_join_todo)->push("TOKEN");
			P->delete_stack(stack_join_todo);
		}

		return;
	}

	current_query_chunk = 0;

	timer.go("Opening the output file");
	Consumer *master_out(options.consumer ? options.consumer : new OutputFile(config.output_file, config.compression == 1));
	if (*output_format == Output_format::daa)
		init_daa(*static_cast<OutputFile*>(master_out));
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
			if (!db_file->load_seqs(&query_block_to_database_id,
				(size_t)(config.chunk_size * 1e9),
				&query_seqs::data_,
				&query_ids::data_,
				true,
				options.db_filter))
				break;
			query_file_offset = db_file->tell_seq();
		}
		else
			if (!load_seqs(query_file->begin(), query_file->end(), *format_n, &query_seqs::data_, query_ids::data_, &query_source_seqs::data_,
				config.store_query_quality ? &query_qual : nullptr,
				(size_t)(config.chunk_size * 1e9), config.qfilt, input_value_traits, paired_mode ? 2 : 1))
				break;

		timer.finish();
		query_seqs::data_->print_stats();

		if (current_query_chunk == 0 && *output_format != Output_format::daa)
			output_format->print_header(*master_out, align_mode.mode, config.matrix.c_str(), score_matrix.gap_open(), score_matrix.gap_extend(), config.max_evalue, query_ids::get()[0],
				unsigned(align_mode.query_translated ? query_source_seqs::get()[0].length() : query_seqs::get()[0].length()));

		if (config.masking == 1 && !options.self) {
			timer.go("Masking queries");
			mask_seqs(*query_seqs::data_, Masking::get());
			timer.finish();
		}

		if (config.multiprocessing)
			P->create_stack_from_file(stack_align_todo, get_ref_part_file_name(stack_align_todo, current_query_chunk));

		run_query_chunk(*db_file, current_query_chunk, *master_out, unaligned_file.get(), aligned_file.get(), metadata, options);

		if (config.multiprocessing)
			P->delete_stack(stack_align_todo);
	}

	if (query_file && !options.query_file) {
		timer.go("Closing the input file");
		query_file->front().close();
		delete query_file;
	}

	timer.go("Closing the output file");
	if (*output_format == Output_format::daa)
		finish_daa(*static_cast<OutputFile*>(master_out), *db_file);
	else
		output_format->print_footer(*master_out);
	master_out->finalize();
	if (!options.consumer) delete master_out;
	if (unaligned_file.get())
		unaligned_file->close();
	if (aligned_file.get())
		aligned_file->close();

	if (!options.db) {
		timer.go("Closing the database file");
		db_file->close();
		delete db_file;
	}

	timer.go("Deallocating taxonomy");
	metadata.free();

	timer.finish();
	log_rss();
	message_stream << "Total time = " << total_timer.get() << "s" << endl;
	statistics.print();
	print_warnings();
}

void run(const Options &options)
{
	task_timer total;

	align_mode = Align_mode(Align_mode::from_command(config.command));

	message_stream << "Temporary directory: " << TempFile::get_temp_dir() << endl;

	if (config.sensitivity >= Sensitivity::VERY_SENSITIVE)
		Config::set_option(config.chunk_size, 0.4);
	else
		Config::set_option(config.chunk_size, 2.0);

	task_timer timer("Opening the database", 1);
	SequenceFile* db_file = options.db ? options.db : SequenceFile::auto_create();
	timer.finish();

	init_output();

	message_stream << "Reference = " << config.database << endl;
	message_stream << "Sequences = " << db_file->sequence_count() << endl;
	message_stream << "Letters = " << db_file->letters() << endl;
	message_stream << "Block size = " << (size_t)(config.chunk_size * 1e9) << endl;
	score_matrix.set_db_letters(config.db_size ? config.db_size : db_file->letters());

	Metadata metadata;
	const bool taxon_filter = !config.taxonlist.empty() || !config.taxon_exclude.empty();
	const bool taxon_culling = config.taxon_k != 0;
	const int db_flags = db_file->metadata();

	int flags = 0;
	if (output_format->needs_taxon_id_lists)
		flags |= SequenceFile::TAXON_MAPPING;
	if (output_format->needs_taxon_nodes)
		flags |= SequenceFile::TAXON_NODES;
	if (output_format->needs_taxon_scientific_names)
		flags |= SequenceFile::TAXON_SCIENTIFIC_NAMES;
	db_file->check_metadata(flags);

	if (output_format->needs_taxon_id_lists || taxon_filter || taxon_culling) {
		if (!(db_flags & SequenceFile::TAXON_MAPPING)) {
			if (taxon_filter)
				throw std::runtime_error("--taxonlist/--taxon-exclude options require taxonomy mapping built into the database.");
			if (taxon_culling)
				throw std::runtime_error("--taxon-k option requires taxonomy mapping built into the database.");
		}
		timer.go("Loading taxonomy mapping");
		metadata.taxon_list = db_file->taxon_list();
		timer.finish();
	}
	if (output_format->needs_taxon_nodes || taxon_filter || taxon_culling) {
		if (!(db_flags & SequenceFile::TAXON_NODES)) {
			if (taxon_filter)
				throw std::runtime_error("--taxonlist/--taxon-exclude options require taxonomy nodes built into the database.");
			if (taxon_culling)
				throw std::runtime_error("--taxon-k option require taxonomy nodes built into the database.");
			if(output_format->needs_taxon_nodes)
				throw std::runtime_error("Output format requires taxonomy nodes built into the database.");
		}
		if (db_file->type() == SequenceFile::DMND && db_file->build_version() < 131) {
			if (taxon_culling)
				throw std::runtime_error("--taxon-k option requires a database built with diamond version >= 0.9.30");
			if (output_format->needs_taxon_ranks)
				throw std::runtime_error("Output fields sskingdoms, skingdoms and sphylums require a database built with diamond version >= 0.9.30");
		}
		timer.go("Loading taxonomy nodes");
		metadata.taxon_nodes = db_file->taxon_nodes();
		if (taxon_filter) {
			timer.go("Building taxonomy filter");
			metadata.taxon_filter = new TaxonomyFilter(config.taxonlist, config.taxon_exclude, *metadata.taxon_list, *metadata.taxon_nodes);
		}
		timer.finish();
	}
	if (output_format->needs_taxon_scientific_names) {
		timer.go("Loading taxonomy names");
		metadata.taxonomy_scientific_names = db_file->taxon_scientific_names();
		timer.finish();
	}

	master_thread(db_file, total, metadata, options);
}

}}
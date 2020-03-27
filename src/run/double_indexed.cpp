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

#include <iostream>
#include <limits>
#include <memory>
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

using namespace std;

namespace Workflow { namespace Search {

static const string reference_partition = "ref_part";

void run_ref_chunk(DatabaseFile &db_file,
	unsigned query_chunk,
	pair<size_t, size_t> query_len_bounds,
	char *query_buffer,
	Consumer &master_out,
	PtrVector<TempFile> &tmp_file,
	const Parameters &params,
	const Metadata &metadata,
	const vector<unsigned> &block_to_database_id)
{
	log_rss();

	task_timer timer;
	if (config.masking == 1) {
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

		for (unsigned i = 0; i < shapes.count(); ++i)
			search_shape(i, query_chunk, query_buffer, ref_buffer);

		timer.go("Deallocating buffers");
		delete[] ref_buffer;
	}

	Consumer* out;
	if (blocked_processing) {
		timer.go("Opening temporary output file");
		tmp_file.push_back(new TempFile());
		out = &tmp_file.back();
	}
	else
		out = &master_out;

	timer.go("Computing alignments");
	align_queries(*Trace_pt_buffer::instance, out, params, metadata);
	delete Trace_pt_buffer::instance;

	if (blocked_processing)
		IntermediateRecord::finish_file(*out);

	timer.go("Deallocating reference");
	delete ref_seqs::data_;
	delete ref_ids::data_;
	timer.finish();
}

void run_query_chunk(DatabaseFile &db_file,
	unsigned query_chunk,
	Consumer &master_out,
	OutputFile *unaligned_file,
	OutputFile *aligned_file,
	const Metadata &metadata,
	const Options &options)
{
	const Parameters params(db_file.ref_header.sequences, db_file.ref_header.letters);
	Parallelizer P = Parallelizer::get();

	task_timer timer("Building query seed set");
	if (query_chunk == 0)
		setup_search_cont();
	if (config.algo == -1) {
		query_seeds = new Seed_set(query_seqs::get(), SINGLE_INDEXED_SEED_SPACE_MAX_COVERAGE);
		timer.finish();
		log_stream << "Seed space coverage = " << query_seeds->coverage() << endl;
		if (use_single_indexed(query_seeds->coverage(), query_seqs::get().letters(), db_file.ref_header.letters))
			config.algo = Config::query_indexed;
		else {
			config.algo = Config::double_indexed;
			delete query_seeds;
			query_seeds = NULL;
		}
	}
	else if (config.algo == Config::query_indexed) {
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
		query_seeds_hashed = new Hashed_seed_set(query_seqs::get());
	}
	timer.finish();

	char *query_buffer = nullptr;
	const pair<size_t, size_t> query_len_bounds = query_seqs::data_->len_bounds(shapes[0].length_);
	setup_search_params(query_len_bounds, 0);

	if (!config.swipe_all) {
		timer.go("Building query histograms");
		query_hst = Partitioned_histogram(*query_seqs::data_, false, &no_filter);

		timer.go("Allocating buffers");
		query_buffer = SeedArray::alloc_buffer(query_hst);
		timer.finish();
	}

	PtrVector<TempFile> tmp_file;
	query_aligned.clear();
	query_aligned.insert(query_aligned.end(), query_ids::get().get_length(), false);
	db_file.rewind();
	vector<unsigned> block_to_database_id;

	if (config.multiprocessing) {
		auto work = P.get_stack(reference_partition);
		string buf;
		while (work->pop(buf)) {
			Chunk chunk = to_chunk(buf);
			// load
			// process
		}
	} else {
		for (current_ref_block = 0;
			 db_file.load_seqs(block_to_database_id, (size_t)(config.chunk_size*1e9),
			 				   &ref_seqs::data_, &ref_ids::data_, true,
								options.db_filter ? options.db_filter : metadata.taxon_filter);
			 ++current_ref_block) {

			run_ref_chunk(db_file, query_chunk, query_len_bounds, query_buffer, master_out, tmp_file, params, metadata, block_to_database_id);
		}
	}

	timer.go("Deallocating buffers");
	delete[] query_buffer;
	delete query_seeds;
	query_seeds = 0;

	log_rss();

	if (blocked_processing) {
		timer.go("Joining output blocks");
		join_blocks(current_ref_block, master_out, tmp_file, params, metadata, db_file);
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
	delete query_seqs::data_;
	delete query_ids::data_;
	delete query_source_seqs::data_;
	delete query_qual;
	if (*output_format != Output_format::daa)
		ReferenceDictionary::get().clear();
}

void master_thread(DatabaseFile *db_file, task_timer &total_timer, Metadata &metadata, const Options &options)
{
	Parallelizer P = Parallelizer::get();
	string mp_reference_partition_file;
	if (config.multiprocessing) {
		P.init();
		mp_reference_partition_file = join_path(P.get_work_directory(), "reference_partition");
		if (P.is_master()) {
			db_file->create_partition((size_t)(config.chunk_size*1e9));
		}
		P.barrier(AUTOTAG);
	}

	task_timer timer("Opening the input file", true);
	TextInputFile *query_file = nullptr;
	const Sequence_file_format *format_n = nullptr;
	if (!options.self) {
		if (config.query_file.empty() && !options.query_file)
			std::cerr << "Query file parameter (--query/-q) is missing. Input will be read from stdin." << endl;
		query_file = options.query_file ? options.query_file : new TextInputFile(config.query_file);
		format_n = guess_format(*query_file);
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

	size_t query_file_offset = 0;

	for (;; ++current_query_chunk) {
		task_timer timer("Loading query sequences", true);

		if (config.multiprocessing) {
			if (P.is_master()) {
				db_file->save_partition(mp_reference_partition_file);
			}
			P.barrier(AUTOTAG);
			P.create_stack_from_file(reference_partition, mp_reference_partition_file);
		}

		if (options.self) {
			db_file->seek_seq(query_file_offset);
			if (!db_file->load_seqs(query_block_to_database_id,
				(size_t)(config.chunk_size * 1e9),
				&query_seqs::data_,
				&query_ids::data_,
				true,
				options.db_filter))
				break;
			query_file_offset = db_file->tell_seq();
		}
		else
			if (!load_seqs(*query_file, *format_n, &query_seqs::data_, query_ids::data_, &query_source_seqs::data_,
				config.store_query_quality ? &query_qual : nullptr,
				(size_t)(config.chunk_size * 1e9), config.qfilt))
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

		run_query_chunk(*db_file, current_query_chunk, *master_out, unaligned_file.get(), aligned_file.get(), metadata, options);
	}

	if (query_file && !options.query_file) {
		timer.go("Closing the input file");
		query_file->close();
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
}

void run(const Options &options)
{
	task_timer total;

	align_mode = Align_mode(Align_mode::from_command(config.command));

	message_stream << "Temporary directory: " << TempFile::get_temp_dir() << endl;

	if (config.mode_very_sensitive) {
		Config::set_option(config.chunk_size, 0.4);
		Config::set_option(config.lowmem, 1u);
	}
	else {
		Config::set_option(config.chunk_size, 2.0);
		Config::set_option(config.lowmem, 4u);
	}

	task_timer timer("Opening the database", 1);
	DatabaseFile *db_file = options.db ? options.db : DatabaseFile::auto_create_from_fasta();
	timer.finish();

	init_output(db_file->has_taxon_id_lists(), db_file->has_taxon_nodes(), db_file->has_taxon_scientific_names());

	message_stream << "Reference = " << config.database << endl;
	message_stream << "Sequences = " << db_file->ref_header.sequences << endl;
	message_stream << "Letters = " << db_file->ref_header.letters << endl;
	message_stream << "Block size = " << (size_t)(config.chunk_size * 1e9) << endl;
	Config::set_option(config.db_size, (uint64_t)db_file->ref_header.letters);
	score_matrix.set_db_letters(db_file->ref_header.letters);

	Metadata metadata;
	const bool taxon_filter = !config.taxonlist.empty() || !config.taxon_exclude.empty();
	const bool taxon_culling = config.taxon_k != 0;
	if (output_format->needs_taxon_id_lists || taxon_filter || taxon_culling) {
		if (db_file->header2.taxon_array_offset == 0) {
			if (taxon_filter)
				throw std::runtime_error("--taxonlist/--taxon-exclude options require taxonomy mapping built into the database.");
			if (taxon_culling)
				throw std::runtime_error("--taxon-k option requires taxonomy mapping built into the database.");
		}
		timer.go("Loading taxonomy mapping");
		metadata.taxon_list = new TaxonList(db_file->seek(db_file->header2.taxon_array_offset), db_file->ref_header.sequences, db_file->header2.taxon_array_size);
		timer.finish();
	}
	if (output_format->needs_taxon_nodes || taxon_filter || taxon_culling) {
		if (db_file->header2.taxon_nodes_offset == 0) {
			if (taxon_filter)
				throw std::runtime_error("--taxonlist/--taxon-exclude options require taxonomy nodes built into the database.");
			if (taxon_culling)
				throw std::runtime_error("--taxon-k option require taxonomy nodes built into the database.");
			if(output_format->needs_taxon_nodes)
				throw std::runtime_error("Output format requires taxonomy nodes built into the database.");
		}
		if (db_file->ref_header.build < 131) {
			if (taxon_culling)
				throw std::runtime_error("--taxon-k option requires a database built with diamond version >= 0.9.30");
			if (output_format->needs_taxon_ranks)
				throw std::runtime_error("Output fields sskingdoms, skingdoms and sphylums require a database built with diamond version >= 0.9.30");
		}
		timer.go("Loading taxonomy nodes");
		metadata.taxon_nodes = new TaxonomyNodes(db_file->seek(db_file->header2.taxon_nodes_offset), db_file->ref_header.build);
		if (taxon_filter) {
			timer.go("Building taxonomy filter");
			metadata.taxon_filter = new TaxonomyFilter(config.taxonlist, config.taxon_exclude, *metadata.taxon_list, *metadata.taxon_nodes);
		}
		timer.finish();
	}
	if (output_format->needs_taxon_scientific_names) {
		timer.go("Loading taxonomy names");
		metadata.taxonomy_scientific_names = new vector<string>;
		db_file->seek(db_file->header2.taxon_names_offset);
		*db_file >> *metadata.taxonomy_scientific_names;
		timer.finish();
	}

	master_thread(db_file, total, metadata, options);
}

}}
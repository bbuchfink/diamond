/****
Copyright (c) 2016, Benjamin Buchfink
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/

#include <numeric>
#include "../data/reference.h"
#include "../output/output_format.h"
#include "../util/seq_file_format.h"
#include "../data/queries.h"
#include "../output/daa_write.h"
#include "../data/load_seqs.h"
#include "../data/index.h"

void search_query_worker(Atomic<unsigned> *next);

void run_query_chunk(Output_stream &master_out)
{
	task_timer timer("Computing alignments");
	Thread_pool threads;
	Atomic<unsigned> query (0);
	for (unsigned i = 0; i < config.threads_; ++i)
		threads.push_back(launch_thread(search_query_worker, &query));
	threads.join_all();
}

void run_ref_chunk(Input_stream &query_file, const Sequence_file_format &input_format, Output_stream &master_out)
{	
	task_timer timer("Building database index");
	shape_from = 0;
	shape_to = 1;
	build_index(ref_seqs::get());
	timer.finish();

	query_file.rewind();

	for (current_query_chunk = 0;; ++current_query_chunk) {
		task_timer timer("Loading query sequences", true);
		size_t n_query_seqs;
		n_query_seqs = load_seqs(query_file, input_format, &query_seqs::data_, query_ids::data_, query_source_seqs::data_, (size_t)(config.chunk_size * 1e9), config.qfilt);
		if (n_query_seqs == 0)
			break;
		timer.finish();
		query_seqs::data_->print_stats();
		run_query_chunk(master_out);
	}

	timer.go("Deallocating memory");
	for (unsigned i = 0; i < shapes.count(); ++i)
		assign_ptr(seed_index[i], new Seed_index());
}

void run_mapper(Database_file &db_file, Timer &total_timer)
{
	task_timer timer("Opening the input file", true);
	auto_ptr<Input_stream> query_file(Compressed_istream::auto_detect(config.query_file));
	const Sequence_file_format *format_n(guess_format(*query_file));

	timer.go("Opening the output file");
	auto_ptr<Output_stream> master_out(config.compression == 1
		? new Compressed_ostream(config.output_file)
		: new Output_stream(config.output_file));
	timer.finish();

	for (current_ref_block = 0; db_file.load_seqs(); ++current_ref_block)
		run_ref_chunk(*query_file, *format_n, *master_out);

	timer.go("Closing the output file");
	master_out->close();
	
	timer.go("Closing the database file");
	db_file.close();

	timer.finish();
	message_stream << "Total wall clock time: " << total_timer.getElapsedTimeInSec() << "s" << endl;
	statistics.print();
}

void run_mapper()
{
	Timer timer2;
	timer2.start();

	Reduction::reduction = Reduction("A KR EDNQ C G H ILVM FYW P ST");

	align_mode = Align_mode(Align_mode::from_command(config.command));
	output_format = auto_ptr<Output_format>(get_output_format());

	message_stream << "Temporary directory: " << Temp_file::get_temp_dir() << endl;

	task_timer timer("Opening the database", 1);
	Database_file db_file;
	timer.finish();
	message_stream << "Reference: " << config.database <<  " (" << ref_header.sequences << " sequences, " << ref_header.letters << " letters)" << endl;
	verbose_stream << "Block size: " << (size_t)(config.chunk_size * 1e9) << endl;
	Config::set_option(config.db_size, (uint64_t)ref_header.letters);

	run_mapper(db_file, timer2);
}

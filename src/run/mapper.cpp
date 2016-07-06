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

void run_mapper(Database_file &db_file, Timer &total_timer)
{
	task_timer timer("Opening the input file", true);
	const Sequence_file_format *format_n(guess_format(config.query_file));
	Compressed_istream query_file(config.query_file);
	current_query_chunk = 0;

	timer.go("Opening the output file");
	auto_ptr<Output_stream> master_out(config.compression == 1
		? new Compressed_ostream(config.output_file)
		: new Output_stream(config.output_file));
	if (*output_format == Output_format::daa)
		init_daa(*master_out);
	timer.finish();

	for (;; ++current_query_chunk) {
		task_timer timer("Loading query sequences", true);
		size_t n_query_seqs;
		n_query_seqs = load_seqs(query_file, *format_n, &query_seqs::data_, query_ids::data_, query_source_seqs::data_, (size_t)(config.chunk_size * 1e9));
		if (n_query_seqs == 0)
			break;
		timer.finish();
		query_seqs::data_->print_stats();

		if (current_query_chunk == 0 && *output_format != Output_format::daa)
			output_format->print_header(*master_out, align_mode.mode, config.matrix.c_str(), config.gap_open, config.gap_extend, config.max_evalue, query_ids::get()[0].c_str(),
				unsigned(align_mode.query_translated ? query_source_seqs::get()[0].length() : query_seqs::get()[0].length()));

		if (align_mode.sequence_type == amino_acid && config.seg == "yes") {
			timer.go("Running complexity filter");
			Complexity_filter::get().run(*query_seqs::data_);
		}

		timer.go("Building query index");
		shape_from = 0;
		shape_to = 1;
		build_query_index();
		
	}

	timer.go("Closing the output file");
	if (*output_format == Output_format::daa)
		finish_daa(*master_out);
	else
		output_format->print_footer(*master_out);
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

	align_mode = Align_mode(Align_mode::from_command(config.command));
	output_format = &get_output_format();

	message_stream << "Temporary directory: " << Temp_file::get_temp_dir() << endl;

	task_timer timer("Opening the database", 1);
	Database_file db_file;
	timer.finish();
	config.set_chunk_size(ref_header.block_size);
	message_stream << "Reference: " << config.database <<  " (" << ref_header.sequences << " sequences, " << ref_header.letters << " letters)" << endl;
	verbose_stream << "Block size: " << (size_t)(ref_header.block_size * 1e9) << endl;
	Config::set_option(config.db_size, (uint64_t)ref_header.letters);

	run_mapper(db_file, timer2);
}

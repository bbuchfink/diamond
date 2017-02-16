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

#include <iostream>
#include <limits>
#include "../data/reference.h"
#include "../data/queries.h"
#include "../basic/statistics.h"
#include "../basic/shape_config.h"
#include "../align/align_queries.h"
#include "../search/align_range.h"
#include "../util/seq_file_format.h"
#include "../data/load_seqs.h"
#include "../output/output_format.h"
#include "../data/frequent_seeds.h"
#include "../output/daa_write.h"

using std::endl;
using std::cout;

struct Search_context
{
	Search_context(unsigned sid, const sorted_list &ref_idx, const sorted_list &query_idx) :
		sid(sid),
		ref_idx(ref_idx),
		query_idx(query_idx)
	{ }
	void operator()(unsigned thread_id, unsigned seedp) const
	{
		Statistics stat;
		align_partition(seedp,
			stat,
			sid,
			ref_idx.get_partition_cbegin(seedp),
			query_idx.get_partition_cbegin(seedp),
			thread_id);
		statistics += stat;
	}
	const unsigned sid;
	const sorted_list &ref_idx;
	const sorted_list &query_idx;
};

void process_shape(unsigned sid,
	Timer &timer_mapping,
	unsigned query_chunk,
	char *query_buffer,
	char *ref_buffer)
{
	using std::vector;

	::partition<unsigned> p(Const::seedp, config.lowmem);
	for (unsigned chunk = 0; chunk < p.parts; ++chunk) {

		message_stream << "Processing query chunk " << query_chunk << ", reference chunk " << current_ref_block << ", shape " << sid << ", index chunk " << chunk << '.' << endl;
		const seedp_range range(p.getMin(chunk), p.getMax(chunk));
		current_range = range;

		task_timer timer("Building reference index", true);
		sorted_list ref_idx(ref_buffer,
			*ref_seqs::data_,
			shapes[sid],
			ref_hst.get(sid),
			range,
			ref_hst.partition());

		timer.go("Building query index");
		timer_mapping.resume();
		sorted_list query_idx(query_buffer,
			*query_seqs::data_,
			shapes[sid],
			query_hst.get(sid),
			range,
			query_hst.partition());

		timer.go("Building seed filter");
		frequent_seeds.build(sid, range, ref_idx, query_idx);

		timer.go("Searching alignments");
		Search_context context(sid, ref_idx, query_idx);
		launch_scheduled_thread_pool(context, Const::seedp, config.threads_);
	}
	timer_mapping.stop();
}

void run_ref_chunk(Database_file &db_file,
	Timer &timer_mapping,
	Timer &total_timer,
	unsigned query_chunk,
	pair<size_t, size_t> query_len_bounds,
	char *query_buffer,
	Output_stream &master_out,
	vector<Temp_file> &tmp_file)
{
	task_timer timer("Building reference histograms");
	const pair<size_t, size_t> len_bounds = ref_seqs::data_->len_bounds(shapes[0].length_);
	ref_hst = Partitioned_histogram(*ref_seqs::data_, (unsigned)len_bounds.second);

	ref_map.init(safe_cast<unsigned>(ref_seqs::get().get_length()));

	timer.go("Allocating buffers");
	char *ref_buffer = sorted_list::alloc_buffer(ref_hst);

	timer.go("Initializing temporary storage");
	timer_mapping.resume();
	Trace_pt_buffer::instance = new Trace_pt_buffer(query_seqs::data_->get_length() / align_mode.query_contexts,
		config.tmpdir,
		config.query_bins);
	timer.finish();
	timer_mapping.stop();

	for (unsigned i = 0; i < shapes.count(); ++i)
		process_shape(i, timer_mapping, query_chunk, query_buffer, ref_buffer);

	timer.go("Deallocating buffers");
	delete[] ref_buffer;

	timer_mapping.resume();
	Output_stream* out;
	if (blocked_processing) {
		timer.go("Opening temporary output file");
		tmp_file.push_back(Temp_file());
		out = new Output_stream(tmp_file.back());
	}
	else
		out = &master_out;

	timer.go("Computing alignments");
	align_queries(*Trace_pt_buffer::instance, out);
	delete Trace_pt_buffer::instance;

	if (blocked_processing) {
		Intermediate_record::finish_file(*out);
		delete out;
	}
	timer_mapping.stop();

	timer.go("Deallocating reference");
	delete ref_seqs::data_;
	delete ref_ids::data_;
	timer.finish();
}

void run_query_chunk(Database_file &db_file,
	Timer &timer_mapping,
	Timer &total_timer,
	unsigned query_chunk,
	pair<size_t, size_t> query_len_bounds,
	Output_stream &master_out,
	Output_stream *unaligned_file)
{
	task_timer timer("Allocating buffers", true);
	char *query_buffer = sorted_list::alloc_buffer(query_hst);
	vector<Temp_file> tmp_file;
	query_aligned.clear();
	query_aligned.insert(query_aligned.end(), query_ids::get().get_length(), false);
	timer.finish();

	db_file.rewind();
	for (current_ref_block = 0; db_file.load_seqs(); ++current_ref_block)
		run_ref_chunk(db_file, timer_mapping, total_timer, query_chunk, query_len_bounds, query_buffer, master_out, tmp_file);

	timer.go("Deallocating buffers");
	timer_mapping.resume();
	delete[] query_buffer;

	if (blocked_processing) {
		timer.go("Joining output blocks");
		join_blocks(current_ref_block, master_out, tmp_file);
	}

	if (unaligned_file) {
		timer.go("Writing unaligned queries");
		write_unaligned(unaligned_file);
	}

	timer.go("Deallocating queries");
	delete query_seqs::data_;
	delete query_ids::data_;
	delete query_source_seqs::data_;
	timer_mapping.stop();
}

void master_thread(Database_file &db_file, Timer &timer_mapping, Timer &total_timer)
{
	task_timer timer("Opening the input file", true);
	timer_mapping.start();
	auto_ptr<Input_stream> query_file(Compressed_istream::auto_detect(config.query_file));
	const Sequence_file_format *format_n(guess_format(*query_file));

	current_query_chunk = 0;

	timer.go("Opening the output file");
	auto_ptr<Output_stream> master_out(config.compression == 1
		? new Compressed_ostream(config.output_file)
		: new Output_stream(config.output_file));
	if (*output_format == Output_format::daa)
		init_daa(*master_out);
	auto_ptr<Output_stream> unaligned_file;
	if (!config.unaligned.empty())
		unaligned_file = auto_ptr<Output_stream>(new Output_stream(config.unaligned));
	timer_mapping.stop();
	timer.finish();

	for (;; ++current_query_chunk) {
		task_timer timer("Loading query sequences", true);
		timer_mapping.resume();
		size_t n_query_seqs;
		n_query_seqs = load_seqs(*query_file, *format_n, &query_seqs::data_, query_ids::data_, query_source_seqs::data_, (size_t)(config.chunk_size * 1e9), config.qfilt);
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

		timer.go("Building query histograms");
		const pair<size_t, size_t> query_len_bounds = query_seqs::data_->len_bounds(shapes[0].length_);
		setup_search_params(query_len_bounds, 0);
		query_hst = Partitioned_histogram(*query_seqs::data_, (unsigned)query_len_bounds.second);
		timer_mapping.stop();
		timer.finish();
		//const bool long_addressing_query = query_seqs::data_->raw_len() > (size_t)std::numeric_limits<uint32_t>::max();

		run_query_chunk(db_file, timer_mapping, total_timer, current_query_chunk, query_len_bounds, *master_out, unaligned_file.get());
	}

	timer.go("Closing the input file");
	query_file->close();

	timer.go("Closing the output file");
	timer_mapping.resume();
	if (*output_format == Output_format::daa)
		finish_daa(*master_out);
	else
		output_format->print_footer(*master_out);
	master_out->close();
	if (unaligned_file.get())
		unaligned_file->close();
	timer_mapping.stop();

	timer.go("Closing the database file");
	db_file.close();

	timer.finish();
	message_stream << "Total time = " << total_timer.getElapsedTimeInSec() << "s" << endl;
	verbose_stream << "Mapping time = " << timer_mapping.getElapsedTimeInSec() << "s" << endl;
	statistics.print();
}

void master_thread_di()
{
	Timer timer2, timer_mapping;
	timer2.start();

	align_mode = Align_mode(Align_mode::from_command(config.command));
	output_format = auto_ptr<Output_format>(get_output_format());

	message_stream << "Temporary directory: " << Temp_file::get_temp_dir() << endl;

	if (config.mode_very_sensitive) {
		Config::set_option(config.chunk_size, 0.4);
		Config::set_option(config.lowmem, 1u);
	}
	else {
		Config::set_option(config.chunk_size, 2.0);
		Config::set_option(config.lowmem, 4u);
	}

	task_timer timer("Opening the database", 1);
	Database_file db_file;
	timer.finish();
	verbose_stream << "Reference = " << config.database << endl;
	verbose_stream << "Sequences = " << ref_header.sequences << endl;
	verbose_stream << "Letters = " << ref_header.letters << endl;
	verbose_stream << "Block size = " << (size_t)(config.chunk_size * 1e9) << endl;
	Config::set_option(config.db_size, (uint64_t)ref_header.letters);

	set_max_open_files(config.query_bins * config.threads_ + unsigned(ref_header.letters / (size_t)(config.chunk_size * 1e9)) + 16);

	master_thread(db_file, timer_mapping, timer2);
}
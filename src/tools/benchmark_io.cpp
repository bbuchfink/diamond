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

#include <string>
#include <random>
#include "../util/system/system.h"
#include "../util/io/output_file.h"
#include "../util/io/input_file.h"
#include "../util/log_stream.h"
#include "../search/trace_pt_buffer.h"
#include "../data/reference.h"
#include "../util/io/input_stream_buffer.h"

using std::vector;
using std::string;
using std::endl;

void benchmark_io() {
	const string file_name = "diamond_io_benchmark.tmp";
	const size_t total_count = 1000000000, query_count = 50;
	
	task_timer timer;

	if (!exists(file_name)) {
		timer.go("Writing output file");
		OutputFile out(file_name);
		std::default_random_engine generator;
		std::uniform_int_distribution<uint32_t> query(0, 2000000), seed(0, 20000), subject(1, UINT32_MAX);
		std::uniform_int_distribution<uint16_t> score(30, 1000);
		for (size_t i = 0; i < total_count / query_count; ++i) {
			out.set(Serializer::VARINT);
			out << query(generator) << seed(generator);
			out.unset(Serializer::VARINT);
			for (size_t j = 0; j < query_count; ++j) {
				out << subject(generator);
				out.write(score(generator));
			}
			out.unset(Serializer::VARINT);
			out << (uint32_t)0;
		}
		const size_t s = out.tell();
		message_stream << "Written " << (double)s / (1 << 30) << "GB. (" << s << ")" << endl;
		message_stream << "Throughput: " << (double)s / (1 << 20) / timer.seconds() << " MB/s" << endl;
		out.close();
	}

	const size_t raw_size = file_size(file_name.c_str());
	message_stream << "File size = " << raw_size << endl;
	timer.go("Reading input file");
	ref_seqs::data_ = new Sequence_set();
	InputFile in(file_name, InputStreamBuffer::ASYNC);
	if (config.raw) {
		char* buf = new char[raw_size];
		in.read(buf, raw_size);
		delete[] buf;
	}
	else {
		vector<hit> out;
		out.reserve(total_count);
		auto it = std::back_inserter(out);
		size_t count = 0;
		try {
			while (true) count += hit::read(in, it);
		}
		catch (EndOfStream&) {}
		message_stream << "Read " << count << " hits." << endl;
	}
	in.close();
	timer.finish();
	delete ref_seqs::data_;	
	message_stream << "Throughput: " << (double)raw_size / (1 << 20) / timer.seconds() << " MB/s" << endl;
}
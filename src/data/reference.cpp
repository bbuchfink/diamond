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

#include <limits>
#include <iostream>
#include <sstream>
#include <set>
#include "../basic/config.h"
#include "reference.h"
#include "../basic/statistics.h"
#include "load_seqs.h"
#include "../util/seq_file_format.h"
#include "../util/log_stream.h"

String_set<0>* ref_ids::data_ = 0;
Ref_map ref_map;
Partitioned_histogram ref_hst;
unsigned current_ref_block;
Reference_header ref_header;
Sequence_set* ref_seqs::data_ = 0;
bool blocked_processing;

using std::cout;
using std::endl;

struct Pos_record
{
	Pos_record()
	{}
	Pos_record(uint64_t pos, size_t len):
		pos(pos),
		seq_len(uint32_t(len))
	{}
	uint64_t pos;
	uint32_t seq_len;
};

void make_db()
{
	message_stream << "Database file: " << config.input_ref_file << endl;
	
	Timer total;
	total.start();
	task_timer timer("Opening the database file", true);
	auto_ptr<Input_stream> db_file (Compressed_istream::auto_detect(config.input_ref_file));
	
	Output_stream out(config.database);
	out.typed_write(&ref_header, 1);

	size_t letters = 0, n = 0;
	uint64_t offset = sizeof(ref_header);

	vector<Letter> seq;
	vector<char> id;
	vector<Pos_record> pos_array;
	FASTA_format format;

	try {

		while (format.get_seq(id, seq, *db_file)) {
			if (seq.size() == 0)
				throw std::runtime_error("File format error: sequence of length 0 at line " + to_string(db_file->line_count));
			if (n % 100000llu == 0llu) {
				std::stringstream ss;
				ss << "Loading sequence data (" << n << " sequences processed)";
				timer.go(ss.str().c_str());
			}
			pos_array.push_back(Pos_record(offset, seq.size()));
			out.write("\xff", 1);
			out.write(seq, false);
			out.write("\xff", 1);
			out.write(id, false);
			out.write("\0", 1);
			letters += seq.size();
			++n;
			offset += seq.size() + id.size() + 3;
		}

	}
	catch (std::exception &e) {
		out.close();
		out.remove();
		throw e;
	}

	timer.go("Writing trailer");
	ref_header.pos_array_offset = offset;
	pos_array.push_back(Pos_record(offset, 0));
	out.write(pos_array, false);

	timer.go("Closing the input file");
	db_file->close();
	
	timer.go("Closing the database file");
	ref_header.letters = letters;
	ref_header.sequences = n;
	out.seekp(0);
	out.typed_write(&ref_header, 1);
	out.close();

	timer.finish();
	message_stream << "Processed " << n << " sequences, " << letters << " letters." << endl;
	message_stream << "Total time = " << total.getElapsedTimeInSec() << "s" << endl;
}

bool Database_file::load_seqs()
{
	task_timer timer("Loading reference sequences");
	const size_t max_letters = (size_t)(config.chunk_size*1e9);
	seek(pos_array_offset);
	size_t letters = 0, seqs = 0, id_letters = 0;

	ref_seqs::data_ = new Sequence_set;
	ref_ids::data_ = new String_set<0>;

	Pos_record r;
	read(&r, 1);
	size_t start_offset = r.pos;

	while (r.seq_len > 0 && letters < max_letters) {
		Pos_record r_next;
		read(&r_next, 1);
		letters += r.seq_len;
		ref_seqs::data_->reserve(r.seq_len);
		const size_t id = r_next.pos - r.pos - r.seq_len - 3;
		id_letters += id;
		ref_ids::data_->reserve(id);
		pos_array_offset += sizeof(Pos_record);
		++seqs;
		r = r_next;
	}

	if (seqs == 0) {
		delete ref_seqs::data_;
		delete ref_ids::data_;
		return false;
	}

	ref_seqs::data_->finish_reserve();
	ref_ids::data_->finish_reserve();
	seek(start_offset);

	for (size_t n = 0; n < seqs; ++n) {
		read(ref_seqs::data_->ptr(n) - 1, ref_seqs::data_->length(n) + 2);
		read(ref_ids::data_->ptr(n), ref_ids::data_->length(n) + 1);
	}

	blocked_processing = seqs < ref_header.sequences;
	return true;
}

void Database_file::get_seq()
{
	vector<Letter> seq;
	string id;
	char c;
	bool all = config.seq_no.size() == 0;
	std::set<size_t> seqs;
	if (!all)
		for (vector<string>::const_iterator i = config.seq_no.begin(); i != config.seq_no.end(); ++i)
			seqs.insert(atoi(i->c_str()) - 1);
	for (size_t n = 0; n < ref_header.sequences; ++n) {
		read(&c, 1);
		while (read(&c, 1), c != '\xff')
			seq.push_back(c);
		while (read(&c, 1), c != '\0')
			id.push_back(c);
		if (all || seqs.find(n) != seqs.end()) {
			cout << '>' << id << endl;
			if (config.reverse) {
				sequence(seq).print(cout, value_traits, sequence::Reversed());
				cout << endl;
			}
			else
				cout << sequence(seq) << endl;
		}
		seq.clear();
		id.clear();
	}
}
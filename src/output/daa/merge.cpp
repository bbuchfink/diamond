/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

#include <unordered_map>
#include "basic/config.h"
#include "daa_file.h"
#include "daa_write.h"
#include "daa_record.h"
#include "util/log_stream.h"

using std::unordered_map;
using std::string;
using std::endl;
using std::vector;

unordered_map<uint32_t, uint32_t> build_mapping(unordered_map<string, uint32_t>& acc2oid, StringSet& seq_ids, vector<uint32_t>& seq_lens, DAAFile& f) {
	TaskTimer timer(("Reading targets for file " + f.file().file_name).c_str());
	unordered_map<uint32_t, uint32_t> r;
	for (uint32_t i = 0; i < f.db_seqs_used(); ++i) {
		auto it = acc2oid.emplace(f.ref_name(i), (uint32_t)acc2oid.size());
		r.emplace(i, it.first->second);
		if (it.second) {
			seq_ids.push_back(f.ref_name(i).begin(), f.ref_name(i).end());
			seq_lens.push_back(f.ref_len(i));
		}
	}
	timer.finish();
	message_stream << "#Targets: " << f.db_seqs_used() << endl;
	return r;
}

static int64_t write_file(DAAFile& f, OutputFile& out, const unordered_map<uint32_t, uint32_t>& subject_map) {
	uint32_t size = 0;
	BinaryBuffer buf;
	TextBuffer out_buf;
	size_t query_num;
	while (f.read_query_buffer(buf, query_num)) {
		DAA_query_record r(f, buf, query_num);
		size_t seek_pos = write_daa_query_record(out_buf, r.query_name.c_str(), r.query_seq.source());
		auto it = r.raw_begin();
		while(it.good()) {
			copy_match_record_raw(it, out_buf, subject_map);
		}
		finish_daa_query_record(out_buf, seek_pos);
		out.write(out_buf.data(), out_buf.size());
		out_buf.clear();
	}
	return int64_t(query_num + 1);
}

void merge_daa() {
	TaskTimer timer("Initializing");
	vector<DAAFile*> files;
	const int n = (int)config.input_ref_file.size();
	if (n == 0)
		throw std::runtime_error("Missing parameter: input files (--in)");
	if (config.output_file.empty())
		throw std::runtime_error("Missing parameter: output file (--out)");
	unordered_map<string, uint32_t> acc2oid;
	vector<unordered_map<uint32_t, uint32_t>> oid_maps;
	StringSet seq_ids;
	vector<uint32_t> seq_lens;
	oid_maps.reserve(n);
	for (int i = 0; i < n; ++i) {
		timer.go("Opening input file");
		files.push_back(new DAAFile(config.input_ref_file[i]));
		timer.finish();
		oid_maps.push_back(build_mapping(acc2oid, seq_ids, seq_lens, *files.back()));
	}
	message_stream << "Total number of targets: " << acc2oid.size() << endl;
	timer.go("Initializing output");
	OutputFile out(config.output_file);
	init_daa(out);
	int64_t query_count = 0;
	for (vector<DAAFile*>::iterator i = files.begin(); i < files.end(); ++i) {
		timer.go(("Writing output for file " + (**i).file().file_name).c_str());
		query_count += write_file(**i, out, oid_maps[i - files.begin()]);
	}
	timer.go("Writing trailer");
	finish_daa(out, *files.front(), seq_ids, seq_lens, query_count);
	out.close();
	timer.finish();
	message_stream << "Total number of queries: " << query_count << endl;
}
/****
DIAMOND protein aligner
Copyright (C) 2022 Max Planck Society for the Advancement of Science e.V.

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

#include <unordered_map>
#include "../../basic/config.h"
#include "daa_file.h"
#include "daa_write.h"

using std::unordered_map;
using std::string;
using std::endl;
using std::vector;

unordered_map<uint32_t, uint32_t> build_mapping(unordered_map<string, uint32_t>& acc2oid, StringSet& seq_ids, vector<uint32_t>& seq_lens, DAA_file& f) {
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

static int64_t write_file(DAA_file& f, OutputFile& out, const unordered_map<uint32_t, uint32_t>& subject_map) {
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
	vector<DAA_file*> files;
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
		files.push_back(new DAA_file(config.input_ref_file[i]));
		timer.finish();
		oid_maps.push_back(build_mapping(acc2oid, seq_ids, seq_lens, *files.back()));
	}
	message_stream << "Total number of targets: " << acc2oid.size() << endl;
	timer.go("Initializing output");
	OutputFile out(config.output_file);
	init_daa(out);
	int64_t query_count = 0;
	for (vector<DAA_file*>::iterator i = files.begin(); i < files.end(); ++i) {
		timer.go(("Writing output for file " + (**i).file().file_name).c_str());
		query_count += write_file(**i, out, oid_maps[i - files.begin()]);
	}
	timer.go("Writing trailer");
	finish_daa(out, *files.front(), seq_ids, seq_lens, query_count);
	out.close();
	timer.finish();
	message_stream << "Total number of queries: " << query_count << endl;
}
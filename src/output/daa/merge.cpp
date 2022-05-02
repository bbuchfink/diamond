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

unordered_map<int64_t, int64_t> build_mapping(unordered_map<string, int64_t>& acc2oid, DAA_file& f) {
	unordered_map<int64_t, int64_t> r;
	for (int64_t i = 0; i < f.db_seqs_used(); ++i) {
		auto it = acc2oid.emplace(f.ref_name(i), acc2oid.size());
		r.emplace(i, it.first->second);
	}
	return r;
}

static void write_file(DAA_file& f, OutputFile& out) {
	uint32_t size = 0;
	BinaryBuffer buf;
	TextBuffer out_buf;
	size_t query_num;
	while (f.read_query_buffer(buf, query_num)) {
		DAA_query_record r(f, buf, query_num);
		size_t seek_pos = write_daa_query_record(out_buf, r.query_name.c_str(), r.query_seq.source());
		DAA_query_record::Match_iterator i = r.begin();
		for (; i.good(); ++i) {
			write_daa_record(out_buf, *i, i->subject_id);
		}
		finish_daa_query_record(out_buf, seek_pos);
		out.write(out_buf.data(), out_buf.size());
		out_buf.clear();
	}
}

void merge_daa() {
	vector<DAA_file*> files;
	const int n = config.input_ref_file.size();
	unordered_map<string, int64_t> acc2oid;
	vector<unordered_map<int64_t, int64_t>> oid_maps;
	oid_maps.reserve(n);
	for (int i = 0; i < n; ++i) {
		files.push_back(new DAA_file(config.input_ref_file[i]));
		oid_maps.push_back(build_mapping(acc2oid, *files.back()));
	}
	OutputFile out(config.output_file);
	init_daa(out);
	for (DAA_file* f : files)
		write_file(*f, out);
	out.close();
}
/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#include "queries.h"

unsigned current_query_chunk;
Sequence_set* query_source_seqs::data_ = 0;
Sequence_set* query_seqs::data_ = 0;
String_set<0>* query_ids::data_ = 0;
Partitioned_histogram query_hst;
vector<bool> query_aligned;
Seed_set *query_seeds = 0;
Hashed_seed_set *query_seeds_hashed = 0;

void write_unaligned(Output_stream *file)
{
	const size_t n = query_ids::get().get_length();
	string s;
	for (size_t i = 0; i < n; ++i) {
		if (!query_aligned[i]) {
			std::stringstream ss;
			ss << '>' << query_ids::get()[i].c_str() << endl;
			(align_mode.query_translated ? query_source_seqs::get()[i] : query_seqs::get()[i]).print(ss, input_value_traits);
			ss << endl;
			s = ss.str();
			file->write(s.data(), s.length());
		}
	}
}
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

#include "global_ranking.h"

namespace Extension { namespace GlobalRanking {

size_t write_merged_query_list_intro(uint32_t query_id, TextBuffer& buf) {
	size_t seek_pos = buf.size();
	buf.write(query_id).write((uint32_t)0);
	return seek_pos;
}

void write_merged_query_list(const IntermediateRecord& r, const ReferenceDictionary& dict, TextBuffer& out, BitVector& ranking_db_filter, Statistics& stat) {
	const uint32_t database_id = dict.database_id(r.subject_dict_id);
	out.write(database_id);
	ranking_db_filter.set(database_id);
	stat.inc(Statistics::TARGET_HITS1);
}

void finish_merged_query_list(TextBuffer& buf, size_t seek_pos) {
	*(uint32_t*)(&buf[seek_pos + sizeof(uint32_t)]) = safe_cast<uint32_t>(buf.size() - seek_pos - sizeof(uint32_t) * 2);
}

}}
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

#ifndef QUERIES_H_
#define QUERIES_H_

#include "../basic/translate.h"
#include "../util/complexity_filter.h"
#include "sorted_list.h"
#include "../basic/statistics.h"
#include "sequence_set.h"
#include "seed_set.h"

extern Partitioned_histogram query_hst;
extern unsigned current_query_chunk;

struct query_source_seqs
{
	static const Sequence_set& get()
	{ return *data_; }
	static Sequence_set *data_;
};

struct query_seqs
{
	static const Sequence_set& get()
	{ return *data_; }
	static Sequence_set *data_;
};

struct query_ids
{
	static const String_set<0>& get()
	{ return *data_; }
	static String_set<0> *data_;
};

extern vector<bool> query_aligned;

void write_unaligned(Output_stream *file);

inline unsigned get_source_query_len(unsigned query_id)
{
	return align_mode.query_translated ? (unsigned)query_seqs::get().reverse_translated_len(query_id*align_mode.query_contexts) : (unsigned)query_seqs::get().length(query_id);
}

extern Seed_set *query_seeds;
extern Hashed_seed_set *query_seeds_hashed;

#endif /* QUERIES_H_ */

/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

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

#ifndef QUERIES_H_
#define QUERIES_H_

#include <mutex>
#include "../basic/translate.h"
#include "../basic/statistics.h"
#include "sequence_set.h"
#include "seed_set.h"
#include "seed_histogram.h"
#include "../util/io/output_file.h"

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
	static Sequence_set& get_nc() {
		return *data_;
	}
	static Sequence_set *data_;
};

struct query_ids
{
	static const String_set<char, 0>& get()
	{ return *data_; }
	static String_set<char, 0> *data_;
};

extern std::mutex query_aligned_mtx;
extern vector<bool> query_aligned;
extern String_set<char, 0> *query_qual;

void write_unaligned(OutputFile *file);
void write_aligned(OutputFile *file);

inline unsigned get_source_query_len(unsigned query_id)
{
	return align_mode.query_translated ? (unsigned)query_seqs::get().reverse_translated_len(query_id*align_mode.query_contexts) : (unsigned)query_seqs::get().length(query_id);
}

inline TranslatedSequence get_translated_query(size_t query_id)
{
	if (align_mode.query_translated)
		return query_seqs::get().translated_seq(query_source_seqs::get()[query_id], query_id*align_mode.query_contexts);
	else
		return TranslatedSequence(query_seqs::get()[query_id]);
}

extern Seed_set *query_seeds;
extern Hashed_seed_set *query_seeds_hashed;
extern vector<unsigned> query_block_to_database_id;

#endif /* QUERIES_H_ */

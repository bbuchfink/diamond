/****
DIAMOND protein aligner
Copyright (C) 2013-2021 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen
						
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

#pragma once

#include "../util/io/text_input_file.h"
#include "sequence_file.h"
#include "seed_histogram.h"

struct ref_seqs
{
	static const SequenceSet& get()
	{
		return *data_;
	}
	static SequenceSet& get_nc()
	{
		return *data_;
	}
	static SequenceSet* data_;
};

struct ref_seqs_unmasked
{
	static const SequenceSet& get()
	{
		return *data_;
	}
	static SequenceSet& get_nc()
	{
		return *data_;
	}
	static SequenceSet* data_;
};

struct ref_ids
{
	static const String_set<char, 0>& get()
	{
		return *data_;
	}
	static String_set<char, 0>* data_;
};

extern Partitioned_histogram ref_hst;
extern unsigned current_ref_block;
extern bool blocked_processing;
extern std::vector<uint32_t> block_to_database_id;

inline size_t max_id_len(const String_set<char, 0>& ids)
{
	size_t max(0);
	for (size_t i = 0; i < ids.get_length(); ++i)
		max = std::max(max, find_first_of(ids[i], Const::id_delimiters));
	return max;
}

inline vector<string> seq_titles(const char* title)
{
	return tokenize(title, "\1");
}

Chunk to_chunk(const string& line);
string to_string(const Chunk& c);

static inline bool long_subject_offsets() {
	return ref_seqs::get().raw_len() > (size_t)std::numeric_limits<uint32_t>::max();
}

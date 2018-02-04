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

#ifndef REFERENCE_H_
#define REFERENCE_H_

#include <memory>
#include <string>
#include <numeric>
#include <limits>
#include "../util/binary_file.h"
#include "sorted_list.h"
#include "../basic/statistics.h"
#include "../data/seed_histogram.h"
#include "../util/hash_table.h"
#include "../util/hash_function.h"
#include "../basic/packed_loc.h"
#include "sequence_set.h"
#include "../util/ptr_vector.h"

using std::auto_ptr;

struct Reference_header
{
	Reference_header() :
		magic_number(0x24af8a415ee186dllu),
		build(Const::build_version),
		db_version(current_db_version),
		sequences(0),
		letters(0)
	{ }
	uint64_t magic_number;
	uint32_t build, db_version;
	uint64_t sequences, letters, pos_array_offset;
	enum { current_db_version = 1 };
};

extern Reference_header ref_header;

struct Database_format_exception : public std::exception
{
	virtual const char* what() const throw()
	{ return "Database file is not a DIAMOND database."; }
};

struct Database_file : public Input_stream
{
	Database_file():
		Input_stream (config.database)
	{
		read_header(*this, ref_header);
		if (ref_header.build < min_build_required || ref_header.db_version != Reference_header::current_db_version)
			throw std::runtime_error("Database was built with a different version of diamond as is incompatible.");
		if (ref_header.sequences == 0)
			throw std::runtime_error("Incomplete database file. Database building did not complete successfully.");
		pos_array_offset = ref_header.pos_array_offset;
	}

	static void read_header(Input_stream &stream, Reference_header &header)
	{
		if (stream.read(&header, 1) != 1)
			throw Database_format_exception();
		if (header.magic_number != Reference_header().magic_number)
			throw Database_format_exception();
	}

	void rewind()
	{
		pos_array_offset = ref_header.pos_array_offset;
	}
	bool load_seqs();
	void get_seq();
	enum { min_build_required = 74 };

	size_t pos_array_offset;
};

void make_db();

struct ref_seqs
{
	static const Sequence_set& get()
	{ return *data_; }
	static Sequence_set& get_nc()
	{ return *data_; }
	static Sequence_set *data_;
};

struct ref_ids
{
	static const String_set<0>& get()
	{ return *data_; }
	static String_set<0> *data_;
};

extern Partitioned_histogram ref_hst;
extern unsigned current_ref_block;
extern bool blocked_processing;

inline size_t max_id_len(const String_set<0> &ids)
{
	size_t max (0);
	for(size_t i=0;i<ids.get_length(); ++i)
		max = std::max(max, find_first_of(ids[i].c_str(), Const::id_delimiters));
	return max;
}

#endif /* REFERENCE_H_ */

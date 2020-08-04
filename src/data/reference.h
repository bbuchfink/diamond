/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
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
#include <vector>
#include <string>
#include <iostream>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include "../util/io/serializer.h"
#include "../util/io/input_file.h"
#include "../util/io/text_input_file.h"
#include "../data/seed_histogram.h"
#include "sequence_set.h"
#include "metadata.h"

using std::vector;
using std::string;

struct ReferenceHeader
{
	ReferenceHeader() :
		magic_number(MAGIC_NUMBER),
		build(Const::build_version),
		db_version(current_db_version),
		sequences(0),
		letters(0)
	{ }
	uint64_t magic_number;
	uint32_t build, db_version;
	uint64_t sequences, letters, pos_array_offset;
	enum { current_db_version = 3 };
	static constexpr uint64_t MAGIC_NUMBER = 0x24af8a415ee186dllu;
	friend InputFile& operator>>(InputFile& file, ReferenceHeader& h);
};

struct ReferenceHeader2
{
	ReferenceHeader2():
		taxon_array_offset(0),
		taxon_array_size(0),
		taxon_nodes_offset(0),
		taxon_names_offset(0)
	{
		memset(hash, 0, sizeof(hash));
	}
	char hash[16];
	uint64_t taxon_array_offset, taxon_array_size, taxon_nodes_offset, taxon_names_offset;

	friend Serializer& operator<<(Serializer &s, const ReferenceHeader2 &h);
	friend Deserializer& operator>>(Deserializer &d, ReferenceHeader2 &h);
};

struct Database_format_exception : public std::exception
{
	virtual const char* what() const throw()
	{ return "Database file is not a DIAMOND database."; }
};


struct Chunk
{
	Chunk() : i(0), offset(0), n_seqs(0)
	{}
	Chunk(size_t i_, size_t offset_, size_t n_seqs_): i(i_), offset(offset_), n_seqs(n_seqs_)
	{}
	size_t i;
	size_t offset;
	size_t n_seqs;
};


struct DatabaseFile : public InputFile
{

	DatabaseFile(const string &file_name);
	DatabaseFile(TempFile &tmp_file);
	static void read_header(InputFile &stream, ReferenceHeader &header);
	static DatabaseFile* auto_create_from_fasta();
	static bool is_diamond_db(const string &file_name);
	void rewind();

	void create_partition(size_t max_letters);
	void create_partition_balanced(size_t max_letters);
	void create_partition_fixednumber(size_t n);

	void save_partition(const string & partition_file_name, const string & annotation = "");
	void load_partition(const string & partition_file_name);
	void clear_partition();
	size_t get_n_partition_chunks();

	bool load_seqs(vector<unsigned> &block_to_database_id, size_t max_letters, Sequence_set **dst_seq, String_set<char, 0> **dst_id, bool load_ids = true, const vector<bool> *filter = NULL, const bool fetch_seqs = true, const Chunk & chunk = Chunk());

	void get_seq();
	void read_seq(string &id, vector<Letter> &seq);
	void skip_seq();
	bool has_taxon_id_lists();
	bool has_taxon_nodes();
	bool has_taxon_scientific_names();
	void close();
	void seek_seq(size_t i);
	size_t tell_seq() const;
	void seek_direct();

	enum { min_build_required = 74, MIN_DB_VERSION = 2 };

	bool temporary;
	size_t pos_array_offset;
	ReferenceHeader ref_header;
	ReferenceHeader2 header2;

	struct Partition
	{
		Partition() : max_letters(0), n_seqs_total(0)
		{}

		size_t max_letters;
		size_t n_seqs_total;
		vector<Chunk> chunks;
	};
	Partition partition;

private:
	void init();

};

void make_db(TempFile **tmp_out = nullptr, TextInputFile *input_file = nullptr);

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
	static const String_set<char, 0>& get()
	{ return *data_; }
	static String_set<char, 0> *data_;
};

extern Partitioned_histogram ref_hst;
extern unsigned current_ref_block;
extern bool blocked_processing;

inline size_t max_id_len(const String_set<char, 0> &ids)
{
	size_t max (0);
	for(size_t i=0;i<ids.get_length(); ++i)
		max = std::max(max, find_first_of(ids[i], Const::id_delimiters));
	return max;
}

inline vector<string> seq_titles(const char *title)
{
	return tokenize(title, "\1");
}

Chunk to_chunk(const string & line);
string to_string(const Chunk & c);

static inline bool long_subject_offsets() {
	return ref_seqs::get().raw_len() > (size_t)std::numeric_limits<uint32_t>::max();
}

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
#include <string>
#include <stdint.h>
#include <limits.h>
#include <list>
#include "../util/io/serializer.h"
#include "../util/io/input_file.h"
#include "../sequence_file.h"

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

struct DatabaseFile : public SequenceFile, public InputFile
{

	DatabaseFile(const string &file_name, Flags flags = Flags::NONE);
	DatabaseFile(TempFile &tmp_file);
	static void read_header(InputFile &stream, ReferenceHeader &header);
	static bool is_diamond_db(const string &file_name);

	void create_partition(size_t max_letters);
	virtual void create_partition_balanced(size_t max_letters) override;
	void create_partition_fixednumber(size_t n);

	virtual void save_partition(const string& partition_file_name, const string& annotation = "") override;
	void load_partition(const string & partition_file_name);
	void clear_partition();
	virtual size_t get_n_partition_chunks() override;
	void skip_seq();
	bool has_taxon_id_lists() const;
	bool has_taxon_nodes() const;
	bool has_taxon_scientific_names() const;
	virtual void close() override;
	virtual void set_seqinfo_ptr(size_t i) override;
	virtual size_t tell_seq() const override;
	virtual void init_seq_access() override;
	static void make_db(TempFile** tmp_out = nullptr, std::list<TextInputFile>* input_file = nullptr);

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

	virtual void init_seqinfo_access() override;
	virtual void seek_chunk(const Chunk& chunk) override;
	virtual SeqInfo read_seqinfo() override;
	virtual void putback_seqinfo() override;
	virtual size_t id_len(const SeqInfo& seq_info, const SeqInfo& seq_info_next) override;
	virtual void seek_offset(size_t p) override;
	virtual void read_seq_data(Letter* dst, size_t len) override;
	virtual void read_id_data(char* dst, size_t len) override;
	virtual void skip_id_data() override;
	virtual size_t sequence_count() const override;
	virtual void read_seq(std::vector<Letter>& seq, std::string& id) override;
	virtual void check_metadata(int flags) const override;
	virtual size_t letters() const override;
	virtual int db_version() const override;
	virtual int program_build_version() const override;
	virtual int metadata() const override;
	virtual TaxonList* taxon_list() override;
	virtual TaxonomyNodes* taxon_nodes() override;
	virtual std::vector<string>* taxon_scientific_names() override;
	virtual int build_version() override;
	virtual ~DatabaseFile();

	static const char* FILE_EXTENSION;

private:
	void init(Flags flags = Flags::NONE);

};

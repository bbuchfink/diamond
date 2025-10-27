/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <string>
#include <stdint.h>
#include "util/io/serializer.h"
#include "util/io/input_file.h"
#include "../sequence_file.h"
#include "../taxon_list.h"
#include "../taxonomy_nodes.h"

struct ReferenceHeader
{
	ReferenceHeader() :
		magic_number(MAGIC_NUMBER),
		build(Const::build_version),
		db_version(current_db_version_prot),
		sequences(0),
		letters(0)
	{ }
	uint64_t magic_number;
	uint32_t build, db_version;
	uint64_t sequences, letters, pos_array_offset;
	static const uint32_t current_db_version_prot;
	static const uint32_t current_db_version_nucl;
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
#ifdef EXTRA
		,db_type(SequenceType::amino_acid)
#endif
	{
		memset(hash, 0, sizeof(hash));
	}
	char hash[16];
	uint64_t taxon_array_offset, taxon_array_size, taxon_nodes_offset, taxon_names_offset;
#ifdef EXTRA
    SequenceType db_type;
#endif

	friend Serializer& operator<<(Serializer &s, const ReferenceHeader2 &h);
	friend Deserializer& operator>>(Deserializer &d, ReferenceHeader2 &h);
};

struct Database_format_exception : public std::exception
{
	virtual const char* what() const throw()
	{ return "Database file is not a DIAMOND database."; }
};

Chunk to_chunk(const std::string& line);
std::string to_string(const Chunk& c);

struct DatabaseFile : public SequenceFile, public InputFile
{

	DatabaseFile(const std::string &file_name, Metadata metadata = Metadata(), Flags flags = Flags::NONE, const ValueTraits& value_traits = amino_acid_traits);
	DatabaseFile(TempFile &tmp_file, const ValueTraits& value_traits = amino_acid_traits);
	static void read_header(InputFile &stream, ReferenceHeader &header);
	static bool is_diamond_db(const std::string &file_name);

	void create_partition(size_t max_letters);
	virtual void create_partition_balanced(int64_t max_letters) override;
	void create_partition_fixednumber(size_t n);

	virtual void save_partition(const std::string& partition_file_name, const std::string& annotation = "") override;
	void load_partition(const std::string & partition_file_name);
	void clear_partition();
	virtual int get_n_partition_chunks() override;
	void skip_seq();
	bool has_taxon_id_lists() const;
	bool has_taxon_nodes() const;
	bool has_taxon_scientific_names() const;
	virtual void close() override;
	virtual void set_seqinfo_ptr(OId i) override;
	virtual OId tell_seq() const override;
	virtual bool eof() const override;
	virtual void init_seq_access() override;
	static void make_db();

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
		std::vector<Chunk> chunks;
	};
	Partition partition;

	virtual int64_t file_count() const override;
	virtual void init_seqinfo_access() override;
	virtual void seek_chunk(const Chunk& chunk) override;
	virtual SeqInfo read_seqinfo() override;
	virtual void putback_seqinfo() override;
	virtual size_t id_len(const SeqInfo& seq_info, const SeqInfo& seq_info_next) override;
	virtual void seek_offset(size_t p) override;
	virtual void read_seq_data(Letter* dst, size_t len, size_t& pos, bool seek) override;
	virtual void read_id_data(const int64_t oid, char* dst, size_t len) override;
	virtual void skip_id_data() override;
	virtual int64_t sequence_count() const override;
	virtual bool read_seq(std::vector<Letter>& seq, std::string& id, std::vector<char>* quals = nullptr) override;
	virtual int64_t letters() const override;
	virtual int db_version() const override;
	virtual int program_build_version() const override;
	virtual std::string taxon_scientific_name(TaxId taxid) const override;
	virtual TaxId get_parent(TaxId taxid) override;
	virtual TaxId max_taxid() const override;
	virtual int rank(TaxId taxid) const override;
	virtual int build_version() override;
	virtual ~DatabaseFile();
	virtual BitVector* filter_by_accession(const std::string& file_name) override;
	virtual const BitVector* builtin_filter() override;
	virtual std::string file_name() override;
	virtual int64_t sparse_sequence_count() const override;
	virtual std::vector<TaxId> taxids(size_t oid) const override;
	virtual void seq_data(size_t oid, std::vector<Letter>& dst) const override;
	virtual size_t seq_length(size_t oid) const override;
	virtual void end_random_access(bool dictionary = true) override;
	virtual void init_write() override;
	virtual void write_seq(const Sequence& seq, const std::string& id) override;

	static const char* FILE_EXTENSION;

private:

	void init(Flags flags = Flags::NONE);
	void read_seqid_list();

	std::unique_ptr<TaxonList> taxon_list_;
	std::vector<std::string> taxon_scientific_names_;
	std::unique_ptr<TaxonomyNodes> taxon_nodes_;

};
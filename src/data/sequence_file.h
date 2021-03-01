/****
DIAMOND protein aligner
Copyright (C) 2021 Max Planck Society for the Advancement of Science e.V.

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

#include "../util/io/input_file.h"
#include "sequence_set.h"
#include "../util/data_structures/bit_vector.h"
#include "taxon_list.h"
#include "taxonomy_nodes.h"
#include "../util/enum.h"

struct Chunk
{
	Chunk() : i(0), offset(0), n_seqs(0)
	{}
	Chunk(size_t i_, size_t offset_, size_t n_seqs_) : i(i_), offset(offset_), n_seqs(n_seqs_)
	{}
	size_t i;
	size_t offset;
	size_t n_seqs;
};

struct SeqInfo
{
	SeqInfo()
	{}
	SeqInfo(uint64_t pos, size_t len) :
		pos(pos),
		seq_len(uint32_t(len))
	{}
	uint64_t pos;
	uint32_t seq_len;
	enum { SIZE = 16 };
};

struct SequenceFile {

	enum class Type { DMND = 0, BLAST = 1 };

	SequenceFile(Type type);

	virtual void init_seqinfo_access() = 0;
	virtual void init_seq_access() = 0;
	virtual void seek_chunk(const Chunk& chunk) = 0;
	virtual size_t tell_seq() const = 0;
	virtual SeqInfo read_seqinfo() = 0;
	virtual void putback_seqinfo() = 0;
	virtual size_t id_len(const SeqInfo& seq_info, const SeqInfo& seq_info_next) = 0;
	virtual void seek_offset(size_t p) = 0;
	virtual void read_seq_data(Letter* dst, size_t len) = 0;
	virtual void read_id_data(char* dst, size_t len) = 0;
	virtual void skip_id_data() = 0;
	virtual size_t sequence_count() const = 0;
	virtual size_t letters() const = 0;
	virtual int db_version() const = 0;
	virtual int program_build_version() const = 0;
	virtual void read_seq(std::vector<Letter>& seq, std::string& id) = 0;
	virtual void check_metadata(int flags) const = 0;
	virtual int metadata() const = 0;
	virtual TaxonList* taxon_list() = 0;
	virtual TaxonomyNodes* taxon_nodes() = 0;
	virtual std::vector<string>* taxon_scientific_names() = 0;
	virtual int build_version() = 0;
	virtual void create_partition_balanced(size_t max_letters) = 0;
	virtual void save_partition(const std::string& partition_file_name, const std::string& annotation = "") = 0;
	virtual size_t get_n_partition_chunks() = 0;
	virtual void set_seqinfo_ptr(size_t i) = 0;
	virtual void close() = 0;

	Type type() const { return type_; }
	bool load_seqs(
		std::vector<uint32_t>* block2db_id,
		const size_t max_letters,
		SequenceSet** dst_seq,
		String_set<char, 0>** dst_id,
		bool load_ids = true,
		const BitVector* filter = nullptr,
		const bool fetch_seqs = true,
		const Chunk& chunk = Chunk());
	void get_seq();
	size_t total_blocks() const;

	static SequenceFile* auto_create(int flags = 0);

	enum Metadata : int {
		TAXON_MAPPING          = 1,
		TAXON_NODES            = 1 << 1,
		TAXON_SCIENTIFIC_NAMES = 1 << 2
	};

	enum {
		NO_COMPATIBILITY_CHECK = 0x1,
		NO_FASTA               = 0x2
	};

private:

	Type type_;

};

template<> struct EnumTraits<SequenceFile::Type> {
	static const EMap<SequenceFile::Type> to_string;
};
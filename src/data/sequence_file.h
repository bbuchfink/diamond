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
#include "../util/data_structures/bit_vector.h"
#include "block.h"

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

struct AccessionNotFound : public std::exception {
};

struct SequenceFile {

	enum class Type { DMND = 0, BLAST = 1 };

	enum class Metadata : int {
		TAXON_MAPPING          = 1,
		TAXON_NODES            = 1 << 1,
		TAXON_SCIENTIFIC_NAMES = 1 << 2,
		TAXON_RANKS            = 1 << 3
	};

	enum class Flags : int {
		NONE                   = 0,
		NO_COMPATIBILITY_CHECK = 1,
		NO_FASTA               = 1 << 1,
		ALL_SEQIDS             = 1 << 2,
		FULL_TITLES            = 1 << 3,
		TARGET_SEQS            = 1 << 4,
		SELF_ALN_SCORES        = 1 << 5
	};

	enum class LoadTitles {
		SINGLE_PASS, LAZY
	};

	SequenceFile(Type type, Alphabet alphabet, Flags flags);

	virtual void init_seqinfo_access() = 0;
	virtual void init_seq_access() = 0;
	virtual void seek_chunk(const Chunk& chunk) = 0;
	virtual size_t tell_seq() const = 0;
	virtual SeqInfo read_seqinfo() = 0;
	virtual void putback_seqinfo() = 0;
	virtual size_t id_len(const SeqInfo& seq_info, const SeqInfo& seq_info_next) = 0;
	virtual void seek_offset(size_t p) = 0;
	virtual void read_seq_data(Letter* dst, size_t len, size_t& pos, bool seek) = 0;
	virtual void read_id_data(const int64_t oid, char* dst, size_t len) = 0;
	virtual void skip_id_data() = 0;
	virtual std::string seqid(size_t oid) const = 0;
	virtual std::string dict_title(size_t dict_id, const size_t ref_block) const = 0;
	virtual size_t dict_len(size_t dict_id, const size_t ref_block) const = 0;
	virtual std::vector<Letter> dict_seq(size_t dict_id, const size_t ref_block) const = 0;
	virtual size_t sequence_count() const = 0;
	virtual size_t sparse_sequence_count() const = 0;
	virtual size_t letters() const = 0;
	virtual int db_version() const = 0;
	virtual int program_build_version() const = 0;
	virtual void read_seq(std::vector<Letter>& seq, std::string& id) = 0;
	virtual Metadata metadata() const = 0;
	virtual TaxonomyNodes* taxon_nodes() = 0;
	virtual std::vector<string>* taxon_scientific_names() = 0;
	virtual int build_version() = 0;
	virtual void create_partition_balanced(size_t max_letters) = 0;
	virtual void save_partition(const std::string& partition_file_name, const std::string& annotation = "") = 0;
	virtual size_t get_n_partition_chunks() = 0;
	virtual void set_seqinfo_ptr(size_t i) = 0;
	virtual void close() = 0;
	virtual void close_weakly() = 0;
	virtual void reopen() = 0;
	virtual BitVector* filter_by_accession(const std::string& file_name) = 0;
	virtual BitVector* filter_by_taxonomy(const std::string& include, const std::string& exclude, TaxonomyNodes& nodes) = 0;
	virtual std::vector<unsigned> taxids(size_t oid) const = 0;
	virtual const BitVector* builtin_filter() = 0;
	virtual std::string file_name() = 0;
	virtual void seq_data(size_t oid, std::vector<Letter>& dst) const = 0;
	virtual size_t seq_length(size_t oid) const = 0;
	virtual void init_random_access(const size_t query_block, const size_t ref_blocks, bool dictionary = true) = 0;
	virtual void end_random_access(bool dictionary = true) = 0;
	virtual std::vector<int> accession_to_oid(const std::string& acc) const = 0;
	virtual LoadTitles load_titles() = 0;
	virtual ~SequenceFile();

	Type type() const { return type_; }
	Block* load_seqs(
		const size_t max_letters,
		bool load_ids = true,
		const BitVector* filter = nullptr,
		bool fetch_seqs = true,
		bool lazy_masking = false,
		const Chunk& chunk = Chunk());
	void get_seq();
	size_t total_blocks() const;
	SequenceSet seqs_by_accession(const std::vector<std::string>::const_iterator begin, const std::vector<std::string>::const_iterator end) const;
	std::vector<Letter> seq_by_accession(const std::string& acc) const;

	void init_dict(const size_t query_block, const size_t target_block);
	void init_dict_block(size_t block, size_t seq_count, bool persist);
	void close_dict_block(bool persist);
	uint32_t dict_id(size_t block, size_t block_id, size_t oid, size_t len, const char* id, const Letter* seq, const double self_aln_score);
	size_t oid(uint32_t dict_id, const size_t ref_block) const;
	double dict_self_aln_score(const size_t dict_size, const size_t ref_block) const;
	size_t dict_size() const {
		return next_dict_id_;
	}
	Flags flags() const {
		return flags_;
	}

	static SequenceFile* auto_create(string& path, Flags flags = Flags::NONE, Metadata metadata = Metadata());

protected:

	void load_dictionary(const size_t query_block, const size_t ref_blocks);
	void free_dictionary();
	static size_t dict_block(const size_t ref_block);

	const Flags flags_;
	std::unique_ptr<OutputFile> dict_file_;
	uint32_t next_dict_id_;
	size_t dict_alloc_size_;
	std::vector<std::vector<uint32_t>> dict_oid_;
	std::vector<std::vector<double>> dict_self_aln_score_;

private:

	enum { DICT_EMPTY = UINT32_MAX };

	virtual void write_dict_entry(size_t block, size_t oid, size_t len, const char* id, const Letter* seq, const double self_aln_score) = 0;
	virtual bool load_dict_entry(InputFile& f, const size_t ref_block) = 0;
	virtual void reserve_dict(const size_t ref_blocks) = 0;
	void load_block(size_t block_id_begin, size_t block_id_end, size_t pos, bool use_filter, const vector<uint64_t>* filtered_pos, bool load_ids, Block* block);
	void load_dict_block(InputFile* f, const size_t ref_block);

	const Type type_;
	const Alphabet alphabet_;

	std::map<size_t, std::vector<uint32_t>> block_to_dict_id_;
	std::mutex dict_mtx_;

};

template<> struct EnumTraits<SequenceFile::Type> {
	static const EMap<SequenceFile::Type> to_string;
};

DEF_ENUM_FLAG_OPERATORS(SequenceFile::Flags)
DEF_ENUM_FLAG_OPERATORS(SequenceFile::Metadata)
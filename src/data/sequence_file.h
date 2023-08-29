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

#include <utility>
#include <unordered_map>
#include "../util/io/input_file.h"
#include "sequence_set.h"
#include "../util/data_structures/bit_vector.h"
#include "taxon_list.h"
#include "taxonomy_nodes.h"
#include "../util/enum.h"
#include "../util/data_structures/bit_vector.h"
#include "block/block.h"
#include "../util/tsv/file.h"

struct Chunk
{
	Chunk() : i(0), offset(0), n_seqs(0)
	{}
	Chunk(int i_, size_t offset_, size_t n_seqs_) : i(i_), offset(offset_), n_seqs(n_seqs_)
	{}
	int i;
	size_t offset;
	int64_t n_seqs;
};

struct AccessionNotFound : public std::exception {
};

struct OperationNotSupported : public std::exception {
};

struct FastaFile;

struct SequenceFile {

	enum class Type { DMND = 0, BLAST = 1, FASTA = 2, BLOCK = 3 };

	enum class Metadata : int {
		TAXON_MAPPING = 1,
		TAXON_NODES = 1 << 1,
		TAXON_SCIENTIFIC_NAMES = 1 << 2,
		TAXON_RANKS = 1 << 3
	};

	enum class Flags : int {
		NONE = 0,
		NO_COMPATIBILITY_CHECK = 1,
		NO_FASTA = 1 << 1,
		ALL_SEQIDS = 1 << 2,
		FULL_TITLES = 1 << 3,
		TARGET_SEQS = 1 << 4,
		SELF_ALN_SCORES = 1 << 5,
		NEED_LETTER_COUNT = 1 << 6,
		ACC_TO_OID_MAPPING = 1 << 7,
		OID_TO_ACC_MAPPING = 1 << 8
	};

	enum class FormatFlags {
		TITLES_LAZY = 1,
		DICT_LENGTHS = 1 << 1,
		DICT_SEQIDS = 1 << 2,
		LENGTH_LOOKUP = 1 << 3,
		SEEKABLE = 1 << 4
	};

	enum class LoadFlags {
		SEQS = 1,
		TITLES = 1 << 1,
		QUALITY = 1 << 2,
		LAZY_MASKING = 1 << 3,
		CONVERT_ALPHABET = 1 << 4,
		NO_CLOSE_WEAKLY = 1 << 5,
        DNA_PRESERVATION = 1 << 6,
        ALL = SEQS | TITLES
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

	SequenceFile(Type type, Alphabet alphabet, Flags flags, FormatFlags format_flags, const ValueTraits& value_traits = amino_acid_traits);

	virtual int64_t file_count() const = 0;
	virtual void init_seqinfo_access() = 0;
	virtual void init_seq_access() = 0;
	virtual void seek_chunk(const Chunk& chunk) = 0;
	virtual OId tell_seq() const = 0;
	virtual bool eof() const = 0;
	virtual bool files_synced();
	virtual SeqInfo read_seqinfo() = 0;
	virtual void putback_seqinfo() = 0;
	virtual size_t id_len(const SeqInfo& seq_info, const SeqInfo& seq_info_next) = 0;
	virtual void seek_offset(size_t p) = 0;
	virtual void read_seq_data(Letter* dst, size_t len, size_t& pos, bool seek) = 0;
	virtual void read_id_data(const int64_t oid, char* dst, size_t len) = 0;
	virtual void skip_id_data() = 0;
	virtual std::string seqid(OId oid) const;
	virtual std::string dict_title(DictId dict_id, const size_t ref_block) const;
	virtual Loc dict_len(DictId dict_id, const size_t ref_block) const;
	virtual std::vector<Letter> dict_seq(DictId dict_id, const size_t ref_block) const;
	virtual int64_t sequence_count() const = 0;
	virtual int64_t sparse_sequence_count() const = 0;
	virtual size_t letters() const = 0;
	virtual int db_version() const = 0;
	virtual int program_build_version() const = 0;
	virtual bool read_seq(std::vector<Letter>& seq, std::string& id, std::vector<char>* quals = nullptr) = 0;
	virtual Metadata metadata() const = 0;
	virtual std::string taxon_scientific_name(TaxId taxid) const;
	virtual int build_version() = 0;
	virtual void create_partition_balanced(size_t max_letters) = 0;
	virtual void save_partition(const std::string& partition_file_name, const std::string& annotation = "") = 0;
	virtual int get_n_partition_chunks() = 0;
	virtual void set_seqinfo_ptr(OId i) = 0;
	virtual void close() = 0;
	virtual void close_weakly() = 0;
	virtual void reopen() = 0;
	virtual BitVector* filter_by_accession(const std::string& file_name) = 0;
	virtual std::vector<TaxId> taxids(size_t oid) const = 0;
	virtual const BitVector* builtin_filter() = 0;
	virtual std::string file_name() = 0;
	virtual void seq_data(size_t oid, std::vector<Letter>& dst) const = 0;
	virtual size_t seq_length(size_t oid) const = 0;
	virtual void init_random_access(const size_t query_block, const size_t ref_blocks, bool dictionary = true) = 0;
	virtual void end_random_access(bool dictionary = true) = 0;
	virtual std::vector<OId> accession_to_oid(const std::string& acc) const;
	virtual void init_write();
	virtual void write_seq(const Sequence& seq, const std::string& id);
	virtual ~SequenceFile();

	Type type() const { return type_; }
	Block* load_seqs(const size_t max_letters, const BitVector* filter = nullptr, LoadFlags flags = LoadFlags(3), const Chunk& chunk = Chunk());
	void get_seq();
	Util::Tsv::File* make_seqid_list();
	size_t total_blocks() const;
	SequenceSet seqs_by_accession(const std::vector<std::string>::const_iterator begin, const std::vector<std::string>::const_iterator end) const;
	std::vector<Letter> seq_by_accession(const std::string& acc) const;
	BitVector* filter_by_taxonomy(const std::string& include, const std::string& exclude) const;
	void write_accession_list(const std::vector<bool>& oids, std::string& file_name);
	template<typename It>
	std::vector<int64_t> seq_offsets(It begin, It end);
	template<typename It>
	FastaFile* sub_db(It begin, It end, const std::string& file_name = std::string());
	template<typename It>
	void sub_db(It begin, It end, FastaFile* out);
	std::vector<std::tuple<FastaFile*, std::vector<OId>, Util::Tsv::File*>> length_sort(int64_t block_size, std::function<int64_t(Loc)>& seq_size);
	Util::Tsv::File& seqid_file();

	void init_dict(const size_t query_block, const size_t target_block);
	void init_dict_block(size_t block, size_t seq_count, bool persist);
	void close_dict_block(bool persist);
	DictId dict_id(size_t block, size_t block_id, size_t oid, size_t len, const char* id, const Letter* seq, const double self_aln_score);
	size_t oid(DictId dict_id, const size_t ref_block) const;
	double dict_self_aln_score(const size_t dict_size, const size_t ref_block) const;
	size_t dict_size() const {
		return next_dict_id_;
	}
	Flags flags() const {
		return flags_;
	}
	FormatFlags format_flags() const {
		return format_flags_;
	}
	const TaxonomyNodes& taxon_nodes() const {
		return *taxon_nodes_;
	}

	static SequenceFile* auto_create(const std::vector<std::string>& path, Flags flags = Flags::NONE, Metadata metadata = Metadata(), const ValueTraits& value_traits = amino_acid_traits);

	size_t mem_size() const {
		size_t n = 0;
		for (auto& v : dict_oid_)
			n += v.size() * sizeof(OId);
		for (auto& v : dict_len_)
			n += v.size() * sizeof(uint32_t);
		for (auto& v : dict_title_)
			n += v.raw_len();
		for (auto& v : dict_seq_)
			n += v.raw_len();
		for (auto& v : dict_self_aln_score_)
			n += v.size() * sizeof(double);
		if (!acc2oid_.empty() || !block_to_dict_id_.empty())
			std::terminate();
		return n;
	}

protected:

	static const char* const SEQID_HDR;

	void load_dictionary(const size_t query_block, const size_t ref_blocks);
	void free_dictionary();
	static size_t dict_block(const size_t ref_block);
	void build_acc_to_oid();
	std::pair<int64_t, int64_t> read_fai_file(const std::string& file_name, int64_t seqs, int64_t letters);
	void add_seqid_mapping(const std::string& id, OId oid);

	const Flags flags_;
	const FormatFlags format_flags_;
	const ValueTraits& value_traits_;
	std::unique_ptr<OutputFile> dict_file_;
	DictId next_dict_id_;
	size_t dict_alloc_size_;
	std::vector<std::vector<OId>> dict_oid_;
	std::vector<std::vector<uint32_t>> dict_len_;
	std::vector<StringSet> dict_title_;
	std::vector<SequenceSet> dict_seq_;
	std::vector<std::vector<double>> dict_self_aln_score_;
	std::unique_ptr<TaxonomyNodes> taxon_nodes_;
	std::unordered_map<std::string, OId> acc2oid_;
	std::unique_ptr<Util::Tsv::File> seqid_file_;
	StringSet acc_;

private:

	static const DictId DICT_EMPTY;

	void write_dict_entry(size_t block, size_t oid, size_t len, const char* id, const Letter* seq, const double self_aln_score);
	bool load_dict_entry(InputFile& f, const size_t ref_block);
	void reserve_dict(const size_t ref_blocks);
	std::pair<Block*, int64_t> load_twopass(const int64_t max_letters, const BitVector* filter, LoadFlags flags, const Chunk& chunk);
	std::pair<Block*, int64_t> load_onepass(const int64_t max_letters, const BitVector* filter, LoadFlags flags);
	void load_dict_block(InputFile* f, const size_t ref_block);

	const Type type_;
	const Alphabet alphabet_;

	std::map<size_t, std::vector<DictId>> block_to_dict_id_;
	std::mutex dict_mtx_;

};

template<> struct EnumTraits<SequenceFile::Type> {
	static const EMap<SequenceFile::Type> to_string;
};

DEFINE_ENUM_FLAG_OPERATORS(SequenceFile::Flags)
DEFINE_ENUM_FLAG_OPERATORS(SequenceFile::Metadata)
DEFINE_ENUM_FLAG_OPERATORS(SequenceFile::FormatFlags)
DEFINE_ENUM_FLAG_OPERATORS(SequenceFile::LoadFlags)

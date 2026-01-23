/****
Copyright © 2013-2026 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

#include <utility>
#include <unordered_map>
#include "util/io/input_file.h"
#include "sequence_set.h"
#include "util/data_structures/bit_vector.h"
#include "util/enum.h"
#include "block/block.h"
#include "util/tsv/tsv.h"

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

constexpr int MAX_LINEAGE = 256;

struct DbFilter {
	DbFilter():
		letter_count(0)
	{ }
	DbFilter(uint64_t size) :
		oid_filter(size),
		letter_count()
	{
	}
	BitVector oid_filter;
	uint64_t letter_count;
};

struct DecodedPackage {
	StringSet ids;
	SequenceSet seqs;
	std::vector<OId> oids;
	std::vector<std::pair<OId, TaxId>> taxids;
	int no;
};

struct RawChunk;

struct SequenceFile {

    enum class Type { DMND = 0, BLAST = 1, FASTA = 2, BLOCK = 3 };

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
		OID_TO_ACC_MAPPING = 1 << 8,
		NEED_LENGTH_LOOKUP = 1 << 9,
		NEED_EARLY_TAXON_MAPPING = 1 << 10,
		TAXON_MAPPING = 1 << 11,
		TAXON_NODES = 1 << 12,
		TAXON_SCIENTIFIC_NAMES = 1 << 13,
		TAXON_RANKS = 1 << 14,
		SEQS = 1 << 15,
		TITLES = 1 << 16,
		QUALITY = 1 << 17,
		LAZY_MASKING = 1 << 18,
		DNA_PRESERVATION = 1 << 19,
		ALL = SEQS | TITLES
	};

	enum class FormatFlags {
		TITLES_LAZY = 1,
		DICT_LENGTHS = 1 << 1,
		DICT_SEQIDS = 1 << 2,
		LENGTH_LOOKUP = 1 << 3,
		SEEKABLE = 1 << 4
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

    SequenceFile(Type type, Flags flags, FormatFlags format_flags, const ValueTraits& value_traits = amino_acid_traits);

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
	virtual void read_id_data(const int64_t oid, char* dst, size_t len, bool all, bool full_titles) = 0;
	virtual void skip_id_data() = 0;
	virtual std::string seqid(OId oid, bool all, bool full_titles);
	virtual std::string dict_title(DictId dict_id, const size_t ref_block) const final;
	virtual Loc dict_len(DictId dict_id, const size_t ref_block);
	virtual std::vector<Letter> dict_seq(DictId dict_id, const size_t ref_block);
	virtual uint64_t sequence_count() const = 0;
	virtual uint64_t letters() const = 0;
	virtual size_t letters_filtered(const DbFilter& v);
	virtual int db_version() const = 0;
	virtual int program_build_version() const = 0;
	virtual bool read_seq(std::vector<Letter>& seq, std::string& id, std::vector<char>* quals = nullptr) = 0;
    Flags metadata() const;
	virtual std::string taxon_scientific_name(TaxId taxid) const;
	virtual int build_version() = 0;
	virtual void create_partition_balanced(int64_t max_letters) = 0;
	virtual void save_partition(const std::string& partition_file_name, const std::string& annotation = "") = 0;
	virtual int get_n_partition_chunks() = 0;
	virtual void set_seqinfo_ptr(OId i) = 0;
	virtual void close() = 0;
	virtual DbFilter* filter_by_accession(const std::string& file_name) = 0;
	virtual std::vector<TaxId> taxids(size_t oid) const = 0;
	virtual TaxId max_taxid() const;
	virtual TaxId get_parent(TaxId taxid);
	virtual int rank(TaxId taxid) const;
	std::set<TaxId> rank_taxid(const std::vector<TaxId>& taxid, int rank);
	TaxId rank_taxid(TaxId taxid, int rank);
	std::vector<TaxId> lineage(TaxId taxid);
	TaxId get_lca(TaxId t1, TaxId t2);
	bool contained(TaxId query, const std::set<TaxId>& filter, bool include_invalid);
	bool contained(const std::vector<TaxId>& query, const std::set<TaxId>& filter, bool all, bool include_invalid);
	virtual std::string file_name() = 0;
	virtual void seq_data(size_t oid, std::vector<Letter>& dst) = 0;
	virtual Loc seq_length(size_t oid) = 0;
	void init_random_access(const size_t query_block, const size_t ref_blocks, bool dictionary = true);
	virtual void end_random_access(bool dictionary = true) = 0;
	virtual std::vector<OId> accession_to_oid(const std::string& acc) const;
	virtual void init_write();
	virtual void write_seq(const Sequence& seq, const std::string& id);
	virtual ~SequenceFile();
	virtual void print_info() const;

	Type type() const { return type_; }
    Block* load_seqs(const int64_t max_letters, const BitVector* filter = nullptr, const Chunk& chunk = Chunk());
	void get_seq();
	Util::Tsv::File* make_seqid_list();
	size_t total_blocks() const;
	SequenceSet seqs_by_accession(const std::vector<std::string>::const_iterator begin, const std::vector<std::string>::const_iterator end);
	std::vector<Letter> seq_by_accession(const std::string& acc);
	DbFilter* filter_by_taxonomy(std::istream& filter, char delimiter, bool exclude);
	void write_accession_list(const std::vector<bool>& oids, std::string& file_name);
	template<typename It>
	std::vector<int64_t> seq_offsets(It begin, It end);
	template<typename It>
	FastaFile* sub_db(It begin, It end, const std::string& file_name = std::string());
	template<typename It>
	void sub_db(It begin, It end, FastaFile* out);
	std::vector<std::tuple<FastaFile*, std::vector<OId>, Util::Tsv::File*>> length_sort(int64_t block_size, std::function<int64_t(Loc)>& seq_size);
	Util::Tsv::File& seqid_file();
	static void init_taxon_output_fields();
    virtual RawChunk* raw_chunk(size_t letters, Flags flags);
	virtual void add_taxid_mapping(const std::vector<std::pair<OId, TaxId>>& taxids);
	virtual int raw_chunk_no() const;

	void init_dict(const size_t query_block, const size_t target_block);
	void init_dict_block(size_t block, size_t seq_count, bool persist);
	void close_dict_block(bool persist);
	DictId dict_id(size_t block, size_t block_id, size_t oid, size_t len, const char* id, const Letter* seq, const double self_aln_score);
	size_t oid(DictId dict_id, const size_t ref_block) const;
	double dict_self_aln_score(const size_t dict_size, const size_t ref_block) const;
	size_t dict_size() const {
		return next_dict_id_;
	}
	Flags& flags() {
		return flags_;
	}
	FormatFlags format_flags() const {
		return format_flags_;
	}

    static SequenceFile* auto_create(const std::vector<std::string>& path, Flags flags = Flags::NONE, const ValueTraits& value_traits = amino_acid_traits);

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
	void init_cache();
	std::pair<Block*, int64_t> load_parallel(const uint64_t max_letters, const BitVector* filter, std::unordered_map<std::string, bool>* accs, const Chunk& chunk, bool load_taxids);

	Flags flags_;
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
	std::unordered_map<std::string, OId> acc2oid_;
	std::unique_ptr<Util::Tsv::File> seqid_file_;
	std::vector<Loc> seq_length_;

private:

	static const DictId DICT_EMPTY;

	void write_dict_entry(size_t block, size_t oid, size_t len, const char* id, const Letter* seq, const double self_aln_score);
	bool load_dict_entry(InputFile& f, const size_t ref_block);
	void reserve_dict(const size_t ref_blocks);
    std::pair<Block*, int64_t> load_twopass(const int64_t max_letters, const BitVector* filter, const Chunk& chunk);
    std::pair<Block*, int64_t> load_onepass(const int64_t max_letters, const BitVector* filter);
	void load_dict_block(InputFile* f, const size_t ref_block);
	void set_cached(TaxId taxon_id, bool contained)
	{
		cached_[taxon_id] = true;
		contained_[taxon_id] = contained;
	}

	const Type type_;

	std::map<size_t, std::vector<DictId>> block_to_dict_id_;
	std::mutex dict_mtx_;
	std::vector<bool> cached_, contained_;

};

template<> struct EnumTraits<SequenceFile::Type> {
	static const EMap<SequenceFile::Type> to_string;
};

struct RawChunk {
	virtual bool empty() const = 0;
	virtual OId begin() const noexcept = 0;
	virtual OId end() const noexcept = 0;
	virtual DecodedPackage* decode(SequenceFile::Flags flags, const BitVector* filter, std::unordered_map<std::string, bool>* accs) const = 0;
	virtual size_t letters() const noexcept = 0;
	virtual size_t bytes() const noexcept = 0;
	virtual ~RawChunk() = default;
	int no;
};

DEFINE_ENUM_FLAG_OPERATORS(SequenceFile::Flags)
DEFINE_ENUM_FLAG_OPERATORS(SequenceFile::FormatFlags)

template<typename It>
void print_taxon_names(It begin, It end, const SequenceFile& db, TextBuffer& out, bool json = false) {
	if (begin == end) {
		out << "N/A";
		return;
	}
	for (It i = begin; i != end; ++i) {
		if (json)
			out << "\"";
		if (i != begin)
			out << (json ? ',' : ';');
		out << db.taxon_scientific_name(*i);
		if (json)
			out << "\"";
	}
}
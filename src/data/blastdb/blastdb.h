#pragma once
#include <objtools/blast/seqdb_reader/seqdbexpert.hpp>
#include <memory>
#include "../sequence_file.h"
#include "../string_set.h"

using BlastOid = int;

struct BlastDB : public SequenceFile {

	BlastDB(const std::string& file_name, Metadata metadata, Flags flags, const ValueTraits& value_traits = amino_acid_traits);

	virtual int64_t file_count() const override;
	virtual void init_seqinfo_access() override;
	virtual void init_seq_access() override;
	virtual void seek_chunk(const Chunk& chunk) override;
	virtual OId tell_seq() const override;
	virtual bool eof() const override;
	virtual SeqInfo read_seqinfo() override;
	virtual void putback_seqinfo() override;
	virtual size_t id_len(const SeqInfo& seq_info, const SeqInfo& seq_info_next) override;
	virtual void seek_offset(size_t p) override;
	virtual void read_seq_data(Letter* dst, size_t len, size_t& pos, bool seek) override;
	virtual void read_id_data(const int64_t oid, char* dst, size_t len) override;
	virtual void skip_id_data() override;
	virtual std::string seqid(OId oid) const override;
	virtual std::string dict_title(DictId dict_id, const size_t ref_block) const override;
	virtual Loc dict_len(DictId dict_id, const size_t ref_block) const override;
	virtual std::vector<Letter> dict_seq(DictId dict_id, const size_t ref_block) const override;
	virtual int64_t sequence_count() const override;
	virtual int64_t sparse_sequence_count() const override;
	virtual size_t letters() const override;
	virtual int db_version() const override;
	virtual int program_build_version() const override;
	virtual bool read_seq(std::vector<Letter>& seq, std::string& id, std::vector<char>* quals = nullptr) override;
	virtual Metadata metadata() const override;
	virtual std::string taxon_scientific_name(TaxId taxid) const override;
	virtual int build_version() override;
	virtual void create_partition_balanced(size_t max_letters) override;
	virtual void save_partition(const std::string& partition_file_name, const std::string& annotation = "") override;
	virtual int get_n_partition_chunks() override;
	virtual void set_seqinfo_ptr(OId i) override;
	virtual void close() override;
	virtual void close_weakly() override;
	virtual void reopen() override;
	virtual BitVector* filter_by_accession(const std::string& file_name) override;
	virtual const BitVector* builtin_filter() override;
	virtual std::string file_name() override;
	virtual std::vector<TaxId> taxids(size_t oid) const override;
	virtual void seq_data(size_t oid, std::vector<Letter>& dst) const override;
	virtual size_t seq_length(size_t oid) const override;
	virtual void init_random_access(const size_t query_block, const size_t ref_blocks, bool dictionary = true) override;
	virtual void end_random_access(bool dictionary = true) override;
	virtual std::vector<OId> accession_to_oid(const std::string& acc) const override;
	virtual void init_write() override;
	virtual void write_seq(const Sequence& seq, const std::string& id) override;
	virtual ~BlastDB();

	static void prep_blast_db(const std::string& path);

	static const char* ACCESSION_FIELD;
	
private:

	const std::string file_name_;
	std::unique_ptr<ncbi::CSeqDBExpert> db_;
	int oid_;
	const bool long_seqids_;
	const Flags flags_;
	int64_t sequence_count_, sparse_sequence_count_;
	BitVector oid_filter_;

	friend void load_blast_seqid();
	friend void load_blast_seqid_lin();

};

#pragma once
#include <objtools/blast/seqdb_reader/seqdbexpert.hpp>
#include <memory>
#include "../sequence_file.h"
#include "../string_set.h"

struct BlastDB : public SequenceFile {

	BlastDB(const std::string& file_name, Metadata metadata, Flags flags);

	virtual void init_seqinfo_access() override;
	virtual void init_seq_access() override;
	virtual void seek_chunk(const Chunk& chunk) override;
	virtual size_t tell_seq() const override;
	virtual SeqInfo read_seqinfo() override;
	virtual void putback_seqinfo() override;
	virtual size_t id_len(const SeqInfo& seq_info, const SeqInfo& seq_info_next) override;
	virtual void seek_offset(size_t p) override;
	virtual void read_seq_data(Letter* dst, size_t len, size_t& pos, bool seek) override;
	virtual void read_id_data(char* dst, size_t len) override;
	virtual void skip_id_data() override;
	virtual std::string seqid(size_t oid) const override;
	virtual std::string dict_title(size_t dict_id, const size_t ref_block) const override;
	virtual size_t dict_len(size_t dict_id, const size_t ref_block) const override;
	virtual std::vector<Letter> dict_seq(size_t dict_id, const size_t ref_block) const override;
	virtual size_t sequence_count() const override;
	virtual size_t sparse_sequence_count() const override;
	virtual size_t letters() const override;
	virtual int db_version() const override;
	virtual int program_build_version() const override;
	virtual void read_seq(std::vector<Letter>& seq, std::string& id) override;
	virtual Metadata metadata() const override;
	virtual TaxonomyNodes* taxon_nodes() override;
	virtual std::vector<string>* taxon_scientific_names() override;
	virtual int build_version() override;
	virtual void create_partition_balanced(size_t max_letters) override;
	virtual void save_partition(const std::string& partition_file_name, const std::string& annotation = "") override;
	virtual size_t get_n_partition_chunks() override;
	virtual void set_seqinfo_ptr(size_t i) override;
	virtual void close() override;
	virtual void close_weakly() override;
	virtual void reopen() override;
	virtual BitVector* filter_by_accession(const std::string& file_name) override;
	virtual BitVector* filter_by_taxonomy(const std::string& include, const std::string& exclude, TaxonomyNodes& nodes) override;
	virtual const BitVector* builtin_filter() override;
	virtual std::string file_name() override;
	virtual std::vector<unsigned> taxids(size_t oid) const override;
	virtual void seq_data(size_t oid, std::vector<Letter>& dst) const override;
	virtual size_t seq_length(size_t oid) const override;
	virtual void init_random_access(const size_t query_block, const size_t ref_blocks, bool dictionary = true) override;
	virtual void end_random_access(bool dictionary = true) override;
	virtual LoadTitles load_titles() override;
	virtual ~BlastDB();

	static const char* ACCESSION_FIELD;
	
private:

	virtual void write_dict_entry(size_t block, size_t oid, size_t len, const char* id, const Letter* seq) override;
	virtual bool load_dict_entry(InputFile& f, const size_t ref_block) override;
	virtual void reserve_dict(const size_t ref_blocks) override;

	const std::string file_name_;
	std::unique_ptr<ncbi::CSeqDBExpert> db_;
	int oid_;
	const bool long_seqids_;
	const Flags flags_;
	BitVector oid_filter_;
	StringSet acc_;

	friend void load_blast_seqid();
	friend void load_blast_seqid_lin();

};
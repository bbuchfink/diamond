#pragma once
#include <objtools/blast/seqdb_reader/seqdbexpert.hpp>
#include "../sequence_file.h"

struct BlastDB : public SequenceFile {

	BlastDB(const std::string& file_name);

	virtual void init_seqinfo_access() override;
	virtual void init_seq_access() override;
	virtual void seek_chunk(const Chunk& chunk) override;
	virtual size_t tell_seq() const override;
	virtual SeqInfo read_seqinfo() override;
	virtual void putback_seqinfo() override;
	virtual size_t id_len(const SeqInfo& seq_info, const SeqInfo& seq_info_next) override;
	virtual void seek_offset(size_t p) override;
	virtual void read_seq_data(Letter* dst, size_t len) override;
	virtual void read_id_data(char* dst, size_t len) override;
	virtual void skip_id_data() override;
	virtual size_t sequence_count() const override;
	virtual size_t letters() const override;
	virtual int db_version() const override;
	virtual int program_build_version() const override;
	virtual void read_seq(std::vector<Letter>& seq, std::string& id) override;
	virtual void check_metadata(int flags) const override;
	virtual int metadata() const override;
	virtual TaxonList* taxon_list() override;
	virtual TaxonomyNodes* taxon_nodes() override;
	virtual std::vector<string>* taxon_scientific_names() override;
	virtual int build_version() override;
	virtual void create_partition_balanced(size_t max_letters) override;
	virtual void save_partition(const std::string& partition_file_name, const std::string& annotation = "") override;
	virtual size_t get_n_partition_chunks() override;
	virtual void set_seqinfo_ptr(size_t i) override;
	virtual void close() override;
	
private:

	ncbi::CSeqDBExpert db_;
	int oid_;
	bool long_seqids_;

};
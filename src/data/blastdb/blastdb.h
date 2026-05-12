/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once
#include "../sequence_file.h"
#include "volume.h"
#include "pal.h"

using BlastOid = int;

struct sqlite3;

struct BlastDB : public SequenceFile {

	BlastDB(const std::string& file_name, Flags flags, const ValueTraits& value_traits = amino_acid_traits);

	virtual void print_info() const override;
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
	virtual void read_id_data(const int64_t oid, char* dst, size_t len, bool all, bool full_titles) override;
	virtual void skip_id_data() override;
	virtual std::string seqid(OId oid, bool all, bool full_titles) override;
	//virtual Loc dict_len(DictId dict_id, const size_t ref_block) override;
	virtual uint64_t sequence_count() const override;
	virtual uint64_t letters() const override;
	virtual int db_version() const override;
	virtual int program_build_version() const override;
	virtual bool read_seq(std::vector<Letter>& seq, std::string& id, std::vector<char>* quals = nullptr) override;
	virtual TaxId max_taxid() const override;
	virtual TaxId get_parent(TaxId taxid) override;
	virtual std::string taxon_scientific_name(TaxId taxid) const override;
	virtual int rank(TaxId taxid) const override {
		auto it = rank_mapping_.find(taxid);
		return it == rank_mapping_.end() ? -1 : it->second;
	}
	virtual int build_version() override;
	virtual void create_partition_balanced(int64_t max_letters) override;
	virtual void save_partition(const std::string& partition_file_name, const std::string& annotation = "") override;
	virtual int get_n_partition_chunks() override;
	virtual void set_seqinfo_ptr(OId i) override;
	virtual void close() override;
	virtual DbFilter* filter_by_accession(const std::string& file_name) override;
	virtual std::string file_name() override;
	virtual std::vector<TaxId> taxids(size_t oid) const override;
	virtual void seq_data(size_t oid, std::vector<Letter>& dst) override;
	virtual Loc seq_length(size_t oid) override;
	virtual void end_random_access(bool dictionary = true) override;
	virtual std::vector<OId> accession_to_oid(const std::string& acc) const override;
	virtual void init_write() override;
	virtual void write_seq(const Sequence& seq, const std::string& id) override;
	virtual ~BlastDB();
	const Pal& pal() const noexcept {
		return pal_;
	}
	virtual RawChunk* raw_chunk(size_t letters, SequenceFile::Flags flags) override;
	virtual void add_taxid_mapping(const std::vector<std::pair<OId, TaxId>>& taxids) override;
	virtual int raw_chunk_no() const override {
		return raw_chunk_no_;
	}

private:

	std::string fetch_seqid(OId oid, bool all, bool fulL_titles);
	std::vector<BlastDefLine> deflines(OId oid, bool all, bool full_titles, bool taxids);
	void open_volume(OId oid);

	const std::string file_name_;
	Pal pal_;
	sqlite3* taxon_db_;
	std::unordered_multimap<OId, TaxId> taxon_mapping_;
	std::map<std::string, int> custom_ranks_;
	std::unordered_map<TaxId, int> rank_mapping_;
	OId oid_;
	const bool long_seqids_;
	Flags flags_;
	BitVector oid_filter_;
	std::vector<TaxId> parent_cache_;
	std::unordered_map<TaxId, std::string> extra_names_;
	std::map<OId, std::string> volumes_;
	BlastVolume volume_;
	int raw_chunk_no_;
	
	friend void load_blast_seqid();
	friend void load_blast_seqid_lin();

};
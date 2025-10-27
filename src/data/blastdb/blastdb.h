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
#include <memory>
#include "../sequence_file.h"

using BlastOid = int;

namespace ncbi {
	class CSeqDB;
	class CLocalTaxon;
}

class sqlite3;

struct BlastDB : public SequenceFile {

	BlastDB(const std::string& file_name, Metadata metadata, Flags flags, const ValueTraits& value_traits = amino_acid_traits);

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
	virtual void read_id_data(const int64_t oid, char* dst, size_t len) override;
	virtual void skip_id_data() override;
	virtual std::string seqid(OId oid, bool all = false) const override;
	virtual Loc dict_len(DictId dict_id, const size_t ref_block) const override;
	virtual std::vector<Letter> dict_seq(DictId dict_id, const size_t ref_block) const override;
	virtual int64_t sequence_count() const override;
	virtual int64_t sparse_sequence_count() const override;
	virtual int64_t letters() const override;
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
	virtual BitVector* filter_by_accession(const std::string& file_name) override;
	virtual const BitVector* builtin_filter() override;
	virtual std::string file_name() override;
	virtual std::vector<TaxId> taxids(size_t oid) const override;
	virtual void seq_data(size_t oid, std::vector<Letter>& dst) const override;
	virtual size_t seq_length(size_t oid) const override;
	virtual void end_random_access(bool dictionary = true) override;
	virtual std::vector<OId> accession_to_oid(const std::string& acc) const override;
	virtual void init_write() override;
	virtual void write_seq(const Sequence& seq, const std::string& id) override;
	virtual ~BlastDB();
	StringSet load_ids(OId begin, OId end) const;

private:

	std::string fetch_seqid(OId oid, bool all) const;

	const std::string file_name_;
	std::unique_ptr<ncbi::CSeqDB> db_;
	//std::unique_ptr<ncbi::CLocalTaxon> taxon_;
	sqlite3* taxon_db_;
	std::map<std::string, int> custom_ranks_;
	std::unordered_map<TaxId, int> rank_mapping_;
	int oid_;
	const bool long_seqids_;
	const Flags flags_;
	int64_t sequence_count_, sparse_sequence_count_;
	BitVector oid_filter_;
	std::vector<TaxId> parent_cache_;
	std::unordered_map<TaxId, std::string> extra_names_;
	std::map<OId, std::string> volumes_;
	
	friend void load_blast_seqid();
	friend void load_blast_seqid_lin();

};
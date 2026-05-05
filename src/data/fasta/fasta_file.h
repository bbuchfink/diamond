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
#include <vector>
#include <list>
#include "../sequence_file.h"
#include "util/tsv/tsv.h"

enum class SeqFileFormat { FASTA, FASTQ };

struct FastaFile : public SequenceFile
{

	struct WriteAccess {};

	FastaFile(const std::vector<std::string> &file_name, Flags flags = Flags::NONE, const ValueTraits& value_traits = amino_acid_traits, const std::string& index_file = std::string());
	FastaFile(const std::string& file_name, bool overwrite, const WriteAccess&, Flags flags = Flags::NONE, const ValueTraits& value_traits = amino_acid_traits);

	virtual int64_t file_count() const override;
	virtual void create_partition_balanced(int64_t max_letters) override;
	virtual void save_partition(const std::string& partition_file_name, const std::string& annotation = "") override;
	virtual int get_n_partition_chunks() override;
	virtual void close() override;
	virtual void set_seqinfo_ptr(OId i) override;
	virtual OId tell_seq() const override;
	virtual bool eof() const override;
	virtual bool files_synced() override;
	virtual void init_seq_access() override;
	virtual void init_seqinfo_access() override;
	virtual void seek_chunk(const Chunk& chunk) override;
	virtual SeqInfo read_seqinfo() override;
	virtual void putback_seqinfo() override;
	virtual size_t id_len(const SeqInfo& seq_info, const SeqInfo& seq_info_next) override;
	virtual void seek_offset(size_t p) override;
	virtual void read_seq_data(Letter* dst, size_t len, size_t& pos, bool seek) override;
	virtual void read_id_data(const int64_t oid, char* dst, size_t len, bool all, bool full_titles) override;
	virtual void skip_id_data() override;
	virtual uint64_t sequence_count() const override;
	virtual bool read_seq(std::vector<Letter>& seq, std::string& id, std::vector<char>* quals = nullptr) override;
	virtual uint64_t letters() const override;
	virtual int db_version() const override;
	virtual int program_build_version() const override;
	virtual int build_version() override;
	virtual ~FastaFile();
	virtual DbFilter* filter_by_accession(const std::string& file_name) override;
	virtual std::string file_name() override;
	virtual std::vector<TaxId> taxids(size_t oid) const override;
	virtual void seq_data(size_t oid, std::vector<Letter>& dst) override;
	virtual Loc seq_length(size_t oid) override;
	virtual void end_random_access(bool dictionary = true) override;
	virtual void init_write() override;
	virtual void write_seq(const Sequence& seq, const std::string& id) override;
	bool is_fasta() const noexcept {
		return format_ == SeqFileFormat::FASTA;
	}
	void advance_seq_count(OId n) {
		oid_ += n;
	}
	virtual RawChunk* raw_chunk(size_t bytes, Flags flags) override;
	virtual int raw_chunk_no() const override {
		return raw_chunk_no_;
	}

private:

	std::pair<int64_t, int64_t> init_read();

	const std::string index_file_;
	std::list<Util::Tsv::File> file_;
	std::list<Util::Tsv::File>::iterator file_ptr_;
	std::unique_ptr<OutputFile> out_file_;
	std::vector<std::streampos> index_;
	SeqFileFormat format_;
	OId oid_;
	int raw_chunk_no_;
	int64_t seqs_, letters_;
	
};
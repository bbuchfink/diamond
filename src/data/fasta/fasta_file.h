/****
Copyright (C) 2012-2026 Benjamin J. Buchfink <buchfink@gmail.com>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

	http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
****/
// SPDX-License-Identifier: Apache-2.0

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
	static void index(const std::string& path, const std::string& dst);

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

private:

	std::pair<int64_t, int64_t> init_read(const std::string& index_file);

	std::list<Util::Tsv::File> file_;
	std::list<Util::Tsv::File>::iterator file_ptr_;
	std::unique_ptr<OutputFile> out_file_;
	SeqFileFormat format_;
	OId oid_;
	int64_t seqs_, letters_;
	
};
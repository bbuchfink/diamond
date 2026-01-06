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
#include <string>
#include "util/ptr_vector.h"
#include "basic/const.h"
#include "util/binary_buffer.h"
#include "util/io/input_file.h"
#include "basic/value.h"

struct DAAFile;

struct DAA_header1
{
	enum { VERSION = 1, COMPATIBILITY_VERSION = 0 };
	DAA_header1():
		magic_number (0x3c0e53476d3ee36bllu),
		version (VERSION)
	{ }
	uint64_t magic_number, version;
};

struct DAA_header2
{
	DAA_header2():
		diamond_build(Const::build_version)
	{
		memset(block_type, 0, sizeof(block_type));
		memset(block_size, 0, sizeof(block_size));
	}
	DAA_header2(uint64_t db_seqs,
			uint64_t db_letters,
			int32_t gap_open,
			int32_t gap_extend,
			int32_t reward,
			int32_t penalty,
			double k,
			double lambda,
			double evalue,
			const std::string &score_matrix,
			unsigned mode):
		diamond_build (Const::build_version),
		db_seqs (db_seqs),
		db_seqs_used (0),
		db_letters (db_letters),
		flags (0),
		query_records (0),
		mode ((int32_t)mode),
		gap_open (gap_open),
		gap_extend (gap_extend),
		reward (reward),
		penalty (penalty),
		reserved1 (0),
		reserved2 (0),
		reserved3 (0),
		k (k),
		lambda (lambda),
		evalue (evalue),
		reserved5 (0)
	{
		memset(block_type, 0, sizeof(block_type));
		memset(block_size, 0, sizeof(block_size));
		memset(this->score_matrix, 0, sizeof(this->score_matrix));
		strcpy(this->score_matrix, score_matrix.c_str());
	}
	DAA_header2(const DAAFile& f);
	typedef enum { empty = 0, alignments = 1, ref_names = 2, ref_lengths = 3 } Block_type;
	uint64_t diamond_build, db_seqs, db_seqs_used, db_letters, flags, query_records;
	int32_t mode, gap_open, gap_extend, reward, penalty, reserved1, reserved2, reserved3;
	double k, lambda, evalue, reserved5;
	char score_matrix[16];
	uint64_t block_size[256];
	char block_type[256];
};

struct DAAFile
{

	DAAFile(const std::string& file_name):
		f_ (file_name),
		query_count_ (0)
	{
		f_.read(&h1_, 1);
		if(h1_.magic_number != DAA_header1().magic_number)
			throw std::runtime_error("Input file is not a DAA file.");
		if(h1_.version > DAA_header1::VERSION)
			throw std::runtime_error("DAA version requires later version of DIAMOND.");
		f_.read(&h2_, 1);

		if(h2_.block_size[0] == 0)
			throw std::runtime_error("Invalid DAA file. DIAMOND run has probably not completed successfully.");

		align_mode = AlignMode(h2_.mode);
		//ref_header.sequences = h2_.db_seqs;

		f_.seek(sizeof(DAA_header1) + sizeof(DAA_header2) + (size_t)h2_.block_size[0]);
		std::string s;
		ref_name_.reserve((size_t)h2_.db_seqs_used);
		for(uint64_t i=0;i<h2_.db_seqs_used;++i) {
			f_ >> s;
			ref_name_.push_back(new std::string(s));
		}
		ref_len_.resize((size_t)h2_.db_seqs_used);
		f_.read(ref_len_.data(), (size_t)h2_.db_seqs_used);

		f_.seek(sizeof(DAA_header1) + sizeof(DAA_header2));
	}

	~DAAFile()
	{
		f_.close();
	}

	uint64_t diamond_build() const
	{ return h2_.diamond_build; }

	uint64_t db_seqs() const
	{ return h2_.db_seqs; }

	uint64_t db_seqs_used() const
	{ return h2_.db_seqs_used; }

	uint64_t db_letters() const
	{ return h2_.db_letters; }

	const char* score_matrix() const
	{ return h2_.score_matrix; }

	int32_t gap_open_penalty() const
	{ return h2_.gap_open; }

	int32_t gap_extension_penalty() const
	{ return h2_.gap_extend; }

	int32_t match_reward() const
	{ return h2_.reward; }

	int32_t mismatch_penalty() const
	{ return h2_.penalty; }

	uint64_t query_records() const
	{ return h2_.query_records; }

	unsigned mode() const
	{ return (unsigned)h2_.mode; }

	const std::string& ref_name(size_t i) const
	{ return ref_name_[i]; }

	uint32_t ref_len(size_t i) const
	{ return ref_len_[i]; }

	double lambda() const
	{
		return h2_.lambda;
	}

	double kappa() const
	{
		return h2_.k;
	}

	double evalue() const
	{
		return h2_.evalue;
	}

	size_t block_size(size_t i) const {
		return h2_.block_size[i];
	}

	const std::vector<uint32_t>& ref_len() const {
		return ref_len_;
	}

	bool read_query_buffer(BinaryBuffer &buf, size_t &query_num)
	{
		uint32_t size;
		f_.read(&size, 1);
		if(size == 0)
			return false;
		buf.clear();
		buf.resize(size);
		f_.read(buf.data(), size);
		query_num = query_count_++;
		return true;
	}

	InputFile& file() {
		return f_;
	}

private:

	InputFile f_;
	size_t query_count_;
	DAA_header1 h1_;
	DAA_header2 h2_;
	PtrVector<std::string> ref_name_;
	std::vector<uint32_t> ref_len_;

	friend void write_file(DAAFile&, OutputFile&);

};
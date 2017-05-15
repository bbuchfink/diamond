/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#ifndef DAA_FILE_H_
#define DAA_FILE_H_

#include <string>
#include <exception>
#include "../util/ptr_vector.h"
#include "../basic/config.h"
#include "../basic/const.h"
#include "../util/binary_buffer.h"
#include "../util/binary_file.h"
#include "../basic/value.h"
#include "../data/reference.h"

using std::string;

struct DAA_header1
{
	DAA_header1():
		magic_number (0x3c0e53476d3ee36bllu),
		version (0)
	{ }
	uint64_t magic_number, version;
};

inline const char* mode_str(int mode)
{
	static const char* mode_str[] = { 0, 0, "blastp", "blastx", "blastn" };
	return mode_str[mode];
}

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
			const string &score_matrix,
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
		strcpy(this->score_matrix, score_matrix.c_str());
	}
	typedef enum { empty = 0, alignments = 1, ref_names = 2, ref_lengths = 3 } Block_type;
	uint64_t diamond_build, db_seqs, db_seqs_used, db_letters, flags, query_records;
	int32_t mode, gap_open, gap_extend, reward, penalty, reserved1, reserved2, reserved3;
	double k, lambda, evalue, reserved5;
	char score_matrix[16];
	uint64_t block_size[256];
	char block_type[256];
};

struct DAA_file
{

	DAA_file(const string& file_name):
		f_ (file_name),
		query_count_ (0)
	{
		f_.read(&h1_, 1);
		if(h1_.magic_number != DAA_header1().magic_number)
			throw std::runtime_error("Input file is not a DAA file.");
		if(h1_.version > Const::daa_version)
			throw std::runtime_error("DAA version requires later version of DIAMOND.");
		f_.read(&h2_, 1);

		if(h2_.block_size[0] == 0)
			throw std::runtime_error("Invalid DAA file. DIAMOND run has probably not completed successfully.");

		Config::set_option(config.db_size, h2_.db_letters);
		align_mode = Align_mode(h2_.mode);
		ref_header.sequences = h2_.db_seqs;

		f_.seek(sizeof(DAA_header1) + sizeof(DAA_header2) + (size_t)h2_.block_size[0]);
		string s;
		ref_name_.reserve((size_t)h2_.db_seqs_used);
		for(uint64_t i=0;i<h2_.db_seqs_used;++i) {
			f_.read_c_str(s);
			ref_name_.push_back(new string(s));
		}
		ref_len_.resize((size_t)h2_.db_seqs_used);
		f_.read(ref_len_.data(), (size_t)h2_.db_seqs_used);

		f_.seek(sizeof(DAA_header1) + sizeof(DAA_header2));
	}

	~DAA_file()
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

	const string& ref_name(size_t i) const
	{ return ref_name_[i]; }

	const uint32_t ref_len(size_t i) const
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

	bool read_query_buffer(Binary_buffer &buf, size_t &query_num)
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

private:

	Input_stream f_;
	size_t query_count_;
	DAA_header1 h1_;
	DAA_header2 h2_;
	Ptr_vector<string> ref_name_;
	vector<uint32_t> ref_len_;

};

#endif /* DAA_FILE_H_ */

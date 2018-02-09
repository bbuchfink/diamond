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

#ifndef SEQ_FILE_FORMAT_H_
#define SEQ_FILE_FORMAT_H_

#include <vector>
#include "../basic/value.h"
#include "io/text_input_file.h"

using std::vector;
using std::pair;

struct file_format_exception : public std::exception
{
	virtual const char* what() const throw()
	{
		return "Invalid input file format";
	}
};

struct Sequence_file_format
{

	virtual bool get_seq(vector<char> &id, vector<Letter> &seq, TextInputFile &s) const = 0;
	virtual ~Sequence_file_format()
	{ }
	
};

struct FASTA_format : public Sequence_file_format
{

	FASTA_format()
	{ }

	virtual bool get_seq(vector<char> &id, vector<Letter> &seq, TextInputFile &s) const;

	virtual ~FASTA_format()
	{ }

};

struct FASTQ_format : public Sequence_file_format
{

	FASTQ_format()
	{ }

	virtual bool get_seq(vector<char> &id, vector<Letter> &seq, TextInputFile &s) const;

	virtual ~FASTQ_format()
	{ }

};

const Sequence_file_format* guess_format(TextInputFile &file);

#endif /* SEQ_FILE_FORMAT_H_ */
/****
Copyright (c) 2014, University of Tuebingen
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****
Author: Benjamin Buchfink
****/

#ifndef SEQ_FILE_FORMAT_H_
#define SEQ_FILE_FORMAT_H_

#include <vector>
#include "../basic/value.h"
#include "compressed_stream.h"

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

	virtual bool get_seq(vector<char> &id, vector<Letter> &seq, Input_stream &s) const = 0;
	virtual ~Sequence_file_format()
	{ }
	
};

struct FASTA_format : public Sequence_file_format
{

	FASTA_format()
	{ }

	virtual bool get_seq(vector<char> &id, vector<Letter> &seq, Input_stream &s) const;

	virtual ~FASTA_format()
	{ }

};

struct FASTQ_format : public Sequence_file_format
{

	FASTQ_format()
	{ }

	virtual bool get_seq(vector<char> &id, vector<Letter> &seq, Input_stream &s) const;

	virtual ~FASTQ_format()
	{ }

};

const Sequence_file_format* guess_format(Input_stream &file);

#endif /* SEQ_FILE_FORMAT_H_ */
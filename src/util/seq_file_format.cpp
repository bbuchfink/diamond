/****
Copyright (c) 2016, Benjamin Buchfink
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/

#include "seq_file_format.h"

struct Raw_text {};
struct Sequence_data {};

template<typename _t>
inline char convert_char(char a)
{
	return a;
}

template<>
inline char convert_char<Sequence_data>(char a)
{
	return input_value_traits.from_char(a);
}

template<typename _t, typename _what>
void copy_line(const string & s, vector<_t>& v, size_t d, _what)
{
	for (string::const_iterator i = s.begin() + d; i != s.end(); ++i)
		v.push_back(convert_char<_what>(*i));
}

bool FASTA_format::get_seq(vector<char>& id, vector<Letter>& seq, Input_stream & s) const
{
	while (s.getline(), s.line.empty() && !s.eof());
	if (s.eof())
		return false;
	if (s.line[0] != '>')
		throw Stream_read_exception(s.line_count, "FASTA format error: Missing '>' at record start.");
	id.clear();
	seq.clear();
	copy_line(s.line, id, 1, Raw_text());
	while (true) {
		s.getline();
		if (s.line.empty()) {
			if (s.eof())
				break;
			else
				continue;
		}
		if (s.line[0] == '>') {
			s.putback_line();
			break;
		}
		try {
			copy_line(s.line, seq, 0, Sequence_data());
		}
		catch (invalid_sequence_char_exception &e) {
			throw Stream_read_exception(s.line_count, e.what());
		}
	}
	return true;
}

bool FASTQ_format::get_seq(vector<char>& id, vector<Letter>& seq, Input_stream & s) const
{
	while (s.getline(), s.line.empty() && !s.eof());
	if (s.eof())
		return false;
	if (s.line[0] != '@')
		throw Stream_read_exception(s.line_count, "FASTQ format error: Missing '@' at record start.");
	id.clear();
	seq.clear();
	copy_line(s.line, id, 1, Raw_text());
	s.getline();
	try {
		copy_line(s.line, seq, 0, Sequence_data());
	}
	catch (invalid_sequence_char_exception &e) {
		throw Stream_read_exception(s.line_count, e.what());
	}
	s.getline();
	if (s.line.empty() || s.line[0] != '+')
		throw Stream_read_exception(s.line_count, "FASTQ format error: Missing '+' line in record.");
	s.getline();
	return true;
}

const Sequence_file_format * guess_format(Input_stream &file)
{
	static const FASTA_format fasta;
	static const FASTQ_format fastq;

	file.getline();
	file.putback_line();
	if (file.line.empty())
		throw std::runtime_error("Error detecting input file format. First line seems to be blank.");
	switch (file.line[0]) {
	case '>': return &fasta;
	case '@': return &fastq;
	default: throw std::runtime_error("Error detecting input file format. First line must begin with '>' (FASTA) or '@' (FASTQ).");
	}
	return 0;
}

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

bool FASTA_format::get_seq(vector<char>& id, vector<Letter>& seq, TextInputFile & s) const
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

bool FASTQ_format::get_seq(vector<char>& id, vector<Letter>& seq, TextInputFile & s) const
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

const Sequence_file_format * guess_format(TextInputFile &file)
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

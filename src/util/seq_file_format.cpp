/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen
						
Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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

#include <memory>
#include "log_stream.h"
#include "seq_file_format.h"
#include "sequence/sequence.h"

using std::unique_ptr;
using std::string;
using std::vector;

struct Raw_text {};
struct Sequence_data {};

template<typename T>
inline char convert_char(char a, const ValueTraits& value_traits)
{
	return a;
}

template<>
inline char convert_char<Sequence_data>(char a, const ValueTraits& value_traits)
{
	return value_traits.from_char(a);
}

template<typename T, typename What>
void copy_line(const string & s, vector<T>& v, size_t d, const ValueTraits& value_traits, What)
{
	for (string::const_iterator i = s.begin() + d; i != s.end(); ++i)
		v.push_back(convert_char<What>(*i, value_traits));
}

bool FASTA_format::get_seq(string& id, vector<Letter>& seq, TextInputFile & s, const ValueTraits& value_traits, vector<char> *qual) const
{
	// !!!
	while (s.getline(), s.line.empty() && !s.eof());
	if (s.line.empty() && s.eof())
		return false;
	if (s.line[0] != '>')
		throw StreamReadException(s.line_count, "FASTA format error: Missing '>' at record start.");
	seq.clear();
	id = s.line.substr(1);
	const char* msg = Util::Seq::fix_title(id);
	if (msg)
		message_stream << "Warning in line " << s.line_count << ": " << msg << std::endl;
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
			copy_line(s.line, seq, 0, value_traits, Sequence_data());
		}
		catch (invalid_sequence_char_exception &e) {
			throw StreamReadException(s.line_count, e.what());
		}
	}
	return true;
}

bool FASTQ_format::get_seq(string& id, vector<Letter>& seq, TextInputFile & s, const ValueTraits& value_traits, vector<char> *qual) const
{
	while (s.getline(), s.line.empty() && !s.eof());
	if (s.line.empty() && s.eof())
		return false;
	if (s.line[0] != '@')
		throw StreamReadException(s.line_count, "FASTQ format error: Missing '@' at record start.");
	seq.clear();
	id = s.line.substr(1);
	s.getline();
	try {
		copy_line(s.line, seq, 0, value_traits, Sequence_data());
	}
	catch (invalid_sequence_char_exception &e) {
		throw StreamReadException(s.line_count, e.what());
	}
	s.getline();
	if (s.line.empty() || s.line[0] != '+')
		throw StreamReadException(s.line_count, "FASTQ format error: Missing '+' line in record.");
	s.getline();
	if (qual) {
		qual->clear();
		copy_line(s.line, *qual, 0, value_traits, Raw_text());
	}
	return true;
}

unique_ptr<const SequenceFileFormat> guess_format(TextInputFile &file)
{
	file.getline();
	file.putback_line();
	if (file.line.empty())
		throw std::runtime_error("Error detecting input file format. First line seems to be blank.");
	switch (file.line[0]) {
	case '>': return unique_ptr<const SequenceFileFormat> { new FASTA_format() };
	case '@': return unique_ptr<const SequenceFileFormat> { new FASTQ_format() };
	default: throw std::runtime_error("Error detecting input file format. First line must begin with '>' (FASTA) or '@' (FASTQ).");
	}
	return 0;
}

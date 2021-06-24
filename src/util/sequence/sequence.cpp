/****
DIAMOND protein aligner
Copyright (C) 2016-2021 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
						
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

#include <algorithm>
#include <stdexcept>
#include "sequence.h"
#include "../util.h"

namespace Util { namespace Seq {

const char* id_delimiters = " \a\b\f\n\r\t\v\1";

void format(Sequence seq, const char *id, const char *qual, OutputFile &out, const std::string &format, const Value_traits &value_traits) {
	static TextBuffer buf;
	constexpr size_t WRAP = 160;
	const size_t l = seq.length();
	if(format == "fasta") {
		buf << '>' << id << '\n';
		for (size_t i = 0; i < seq.length(); i += WRAP) {
			seq.print(buf, i, std::min(i + WRAP, seq.length()), value_traits);
			buf << '\n';
		}
	}
	else if (format == "fastq") {
		buf << '@' << id << '\n';
		seq.print(buf, value_traits);
		buf << "\n+\n";
		buf << qual << '\n';
	}
	else
		throw std::runtime_error("Invalid sequence file format");
	out.write(buf.data(), buf.size());
	buf.clear();
}

std::string all_seqids(const char* s)
{
	string r;
	const vector<string> t(tokenize(s, "\1"));
	for (vector<string>::const_iterator i = t.begin(); i != t.end(); ++i) {
		if (i != t.begin())
			r.append("\1");
		r.append(i->substr(0, find_first_of(i->c_str(), id_delimiters)));
	}
	return r;
}

std::string seqid(const char* title)
{
	return string(title, title + find_first_of(title, id_delimiters));
}

void get_title_def(const std::string& s, std::string& title, std::string& def)
{
	const size_t i = find_first_of(s.c_str(), id_delimiters);
	title = s.substr(0, i);
	if (i >= s.length())
		def.clear();
	else
		def = s.substr(i + 1);
}

}}
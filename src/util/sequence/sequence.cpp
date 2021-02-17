/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

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

using namespace std;

namespace Util { namespace Seq {

void format(Sequence seq, const char *id, const char *qual, OutputFile &out, const std::string &format, const Value_traits &value_traits) {
	static TextBuffer buf;
	constexpr size_t WRAP = 160;
	const size_t l = seq.length();
	if(format == "fasta") {
		buf << '>' << id << '\n';
		for (size_t i = 0; i < seq.length(); i += WRAP) {
			seq.print(buf, i, min(i + WRAP, seq.length()), value_traits);
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
		throw runtime_error("Invalid sequence file format");
	out.write(buf.get_begin(), buf.size());
	buf.clear();
}

}}
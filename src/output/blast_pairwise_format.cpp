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

#include "output_format.h"

void Pairwise_format::print_match(const Hsp_context& r, Text_buffer &out)
{
	static const unsigned width = 60;
	out << '>';
	Output_format::print_title(out, r.subject_name, true, true, " ");
	out << "\nLength=" << r.subject_len << "\n\n";
	out << " Score = " << r.bit_score() << " bits (" << r.score() << "),  Expect = ";
	out.print_e(r.evalue());
	out << '\n';
	out << " Identities = " << r.identities() << '/' << r.length() << " (" << percentage<unsigned,unsigned>(r.identities(), r.length()) << "%), Positives = " << r.positives() << '/' << r.length() << " (" << percentage<unsigned, unsigned>(r.positives(), r.length())
		<< "%), Gaps = " << r.gaps() << '/' << r.length() << " (" << percentage<unsigned, unsigned>(r.gaps(), r.length()) << "%)\n";
	if (align_mode.query_translated)
		out << " Frame = " << r.blast_query_frame() << '\n';
	out << '\n';
	const unsigned digits = (unsigned)std::max(ceil(log10(r.subject_range().end_)), ceil(log10(r.query_range().end_)));

	Hsp_context::Iterator qi = r.begin(), mi = r.begin(), si = r.begin();
	while (qi.good()) {
		out << "Query  ";
		out.print(qi.query_pos+1, digits);
		out << "  ";
		for (unsigned i = 0; i < width && qi.good(); ++i, ++qi)
			out << qi.query_char();
		out << " ";
		out.print(qi.query_pos, 0);
		out << '\n';

		for (unsigned i = 0; i < digits + 9; ++i)
			out << ' ';
		for (unsigned i = 0; i < width && mi.good(); ++i, ++mi)
			out << mi.midline_char();
		out << '\n';

		out << "Sbjct  ";
		out.print(si.subject_pos+1, digits);
		out << "  ";
		for (unsigned i = 0; i < width && si.good(); ++i, ++si)
			out << si.subject_char();
		out << " ";
		out.print(si.subject_pos, 0);
		out << "\n\n";
	}
}

void Pairwise_format::print_footer(Output_stream &out) const
{

}

void Pairwise_format::print_query_epilog(Text_buffer &out, const char *query_title, bool unaligned) const
{
}

void Pairwise_format::print_query_intro(size_t query_num, const char *query_name, unsigned query_len, Text_buffer &out, bool unaligned) const
{
	out << "Query= " << query_name << "\n\nLength=" << query_len << "\n\n";
	if (unaligned) {
		out << "\n***** No hits found *****\n\n\n";
	}
}

void Pairwise_format::print_header(Output_stream &f, int mode, const char *matrix, int gap_open, int gap_extend, double evalue, const char *first_query_name, unsigned first_query_len) const
{
	static const char* header = "BLASTP 2.3.0+\n\n\n";
	f.write(header, strlen(header));
}
/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

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

#include "output_format.h"
#include "../util/util.h"

void Pairwise_format::print_match(const HspContext& r, Output::Info &info)
{
	static const unsigned width = 60;
	TextBuffer& out = info.out;
	const int dna_len = (int)r.query.source().length();
	const Strand strand = r.frame() < 3 ? FORWARD : REVERSE;
	out << '>';
	OutputFormat::print_title(out, r.target_title.c_str(), true, true, " ");
	out << "\nLength=" << r.subject_len << "\n\n";
	out << " Score = " << r.bit_score() << " bits (" << r.score() << "),  Expect = ";
	out.print_e(r.evalue());
	out << '\n';
	out << " Identities = " << r.identities() << '/' << r.length() << " (" << percentage<unsigned,unsigned>(r.identities(), r.length()) << "%), Positives = " << r.positives() << '/' << r.length() << " (" << percentage<unsigned, unsigned>(r.positives(), r.length())
		<< "%), Gaps = " << r.gaps() << '/' << r.length() << " (" << percentage<unsigned, unsigned>(r.gaps(), r.length()) << "%)\n";
	if (align_mode.query_translated)
		out << " Frame = " << r.blast_query_frame() << '\n';
	out << '\n';
	const unsigned digits = (unsigned)std::max(ceil(log10(r.subject_range().end_)), ceil(log10(r.query_source_range().end_)));

	HspContext::Iterator qi = r.begin(), mi = r.begin(), si = r.begin();
	while (qi.good()) {
		out << "Query  ";
		out.print(qi.query_pos.absolute(dna_len) + 1, digits);
		out << "  ";
		for (unsigned i = 0; i < width && qi.good(); ++i, ++qi)
			out << qi.query_char();
		out << " ";
		out << TranslatedPosition::oriented_position(qi.query_pos.in_strand() - 1, strand, dna_len) + 1;
		out << '\n';

		for (unsigned i = 0; i < digits + 9; ++i)
			out << ' ';
		for (unsigned i = 0; i < width && mi.good(); ++i, ++mi)
			out << mi.midline_char(score_matrix(mi.query(), mi.subject()));
		out << '\n';

		out << "Sbjct  ";
		out.print(si.subject_pos+1, digits);
		out << "  ";
		for (unsigned i = 0; i < width && si.good(); ++i, ++si)
			out << si.subject_char();
		out << " ";
		out << si.subject_pos;
		out << "\n\n";
	}
}

void Pairwise_format::print_footer(Consumer &out) const
{

}

void Pairwise_format::print_query_epilog(Output::Info &info) const
{
}

void Pairwise_format::print_query_intro(Output::Info &info) const
{
	info.out << "Query= " << info.query.title << "\n\nLength=" << info.query.len << "\n\n";
	if (info.unaligned) {
		info.out << "\n***** No hits found *****\n\n\n";
	}
}

void Pairwise_format::print_header(Consumer &f, int mode, const char *matrix, int gap_open, int gap_extend, double evalue, const char *first_query_name, unsigned first_query_len) const
{
	static const char* header = "BLASTP 2.3.0+\n\n\n";
	f.consume(header, strlen(header));
}
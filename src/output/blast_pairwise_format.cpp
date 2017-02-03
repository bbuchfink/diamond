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

#include "output_format.h"

void Pairwise_format::print_match(const Hsp_context& r, Text_buffer &out) const
{
	static const unsigned width = 60;
	out << '>' << r.subject_name << '\n';
	out << "Length=" << r.subject_len << "\n\n";
	out << " Score = " << r.bit_score() << " bits (" << r.score() << "),  Expect = ";
	out.print_e(r.evalue());
	out << '\n';
	out << " Identities = " << r.identities() << '/' << r.length() << " (" << percentage<unsigned,unsigned>(r.identities(), r.length()) << "%), Positives = " << r.positives() << '/' << r.length() << " (" << percentage<unsigned, unsigned>(r.positives(), r.length())
		<< "%), Gaps = " << r.gaps() << '/' << r.length() << " (" << percentage<unsigned, unsigned>(r.gaps(), r.length()) << "%)\n";
	if (align_mode.query_translated)
		out << " Frame = " << r.blast_query_frame() << '\n';
	out << '\n';

	Hsp_context::Iterator qi = r.begin(), mi = r.begin(), si = r.begin();
	while (qi.good()) {
		out << "Query  ";
		out.print(qi.query_pos+1, 0);
		out << "  ";
		for (unsigned i = 0; i < width && qi.good(); ++i, ++qi)
			out << qi.query_char();
		out << " ";
		out.print(qi.query_pos, 0);
		out << '\n';

		out << "             ";
		for (unsigned i = 0; i < width && mi.good(); ++i, ++mi)
			out << mi.midline_char();
		out << '\n';

		out << "Sbjct  ";
		out.print(si.subject_pos+1, 0);
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

void Pairwise_format::print_query_epilog(Text_buffer &out, bool unaligned) const
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
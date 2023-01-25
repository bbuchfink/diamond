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

#include <sstream>
#include "output_format.h"
#include "../util/sequence/sequence.h"

using namespace std;

void print_md(const HspContext &r, TextBuffer &buf)
{
	unsigned matches = 0, del = 0;
	for (Packed_transcript::Const_iterator i = r.begin_old(); i.good(); ++i) {
		switch (i->op) {
		case op_match:
			del = 0;
			matches += i->count;
			break;
		case op_substitution:
			if (matches > 0) {
				buf << matches;
				matches = 0;
			}
			else if (del > 0) {
				buf << '0';
				del = 0;
			}
			buf << value_traits.alphabet[(long)i->letter];
			break;
		case op_deletion:
			if (matches > 0) {
				buf << matches;
				matches = 0;
			}
			if (del == 0)
				buf << '^';
			buf << value_traits.alphabet[(long)i->letter];
			++del;
		default:
			break;
		}
	}
	if (matches > 0)
		buf << matches;
}

void print_cigar(const HspContext &r, TextBuffer &buf)
{
	static const unsigned map[] = { 0, 1, 2, 0, 3, 4 };
	static const char letter[] = { 'M', 'I', 'D', '\\', '/' };
	unsigned n = 0, op = 0;
	for (Packed_transcript::Const_iterator i = r.begin_old(); i.good(); ++i) {
		if (map[i->op] == op)
			n += i->count;
		else {
			if (n > 0)
				buf << n << letter[op];
			n = i->count;
			op = map[i->op];
		}
	}
	if (n > 0)
		buf << n << letter[op];
}

void Sam_format::print_query_intro(Output::Info &info) const
{
	if (info.unaligned) {
		info.out.write_until(info.query.title, Util::Seq::id_delimiters);
		info.out << "\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n";
	}
}

void Sam_format::print_match(const HspContext& r, Output::Info &info)
{
	TextBuffer& out = info.out;
	out.write_until(r.query_title.c_str(), Util::Seq::id_delimiters);
	out << '\t' << '0' << '\t';

	const bool lt = (config.salltitles || (config.command == Config::view)) ? true : false;
	print_title(out, r.target_title.c_str(), lt, lt, "<>");

	out << '\t'
		<< r.subject_range().begin_ + 1 << '\t'
		<< "255" << '\t';

	print_cigar(r, out);

	out << '\t'
		<< '*' << '\t'
		<< '0' << '\t'
		<< '0' << '\t'
		<< Sequence(r.query.index(r.frame()).data() + r.query_range().begin_, r.query_range().length()) << '\t'
		<< '*' << '\t'
		<< "AS:i:" << (uint32_t)score_matrix.bitscore(r.score()) << '\t'
		<< "NM:i:" << r.length() - r.identities() << '\t'
		<< "ZL:i:" << r.subject_len << '\t'
		<< "ZR:i:" << r.score() << '\t'
		<< "ZE:f:";
	out.print_e(r.evalue());
	out << '\t'
		<< "ZI:i:" << r.identities() * 100 / r.length() << '\t'
		<< "ZF:i:" << Frame(r.frame()).signed_frame() << '\t'
		<< "ZS:i:" << r.oriented_query_range().begin_ + 1 << '\t'
		<< "MD:Z:";

	print_md(r, out);
	out << '\n';
}

void Sam_format::print_header(Consumer &f, int mode, const char *matrix, int gap_open, int gap_extend, double evalue, const char *first_query_name, unsigned first_query_len) const
{
	static const char* mode_str[] = { 0, 0, "BlastP", "BlastX", "BlastN" };
	stringstream h;
	h << "@HD\tVN:1.5\tSO:query" << endl;
	h << "@PG\tPN:DIAMOND\tVN:" << Const::version_string << "\tCL:" << config.invocation << endl;
	h << "@mm\t" << mode_str[mode] << endl;
	h << "@CO\t" << mode_str[mode] << "-like alignments" << endl;
	h << "@CO\tReporting AS: bitScore, ZR: rawScore, ZE: expected, ZI: percent identity, ZL: reference length, ZF: frame, ZS: query start DNA coordinate" << endl;

	const string s(h.str());
	f.consume(s.c_str(), s.length());
}
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

#ifndef OUTPUT_FORMAT_H_
#define OUTPUT_FORMAT_H_

#include <exception>
#include "../basic/match.h"
#include "../align/match_func.h"
#include "../output/daa_file.h"
#include "../output/daa_record.h"

struct Output_format
{
	virtual void print_match(const DAA_match_record &r, Text_buffer &out) const = 0;
	virtual void print_header(Compressed_ostream &f, int mode) const
	{ }
	virtual ~Output_format()
	{ }
	static size_t print_salltitles(Text_buffer &buf, const char *id)
	{
		size_t n = 0;
		const vector<string> t (tokenize(id, "\1"));
		vector<string>::const_iterator i=t.begin();
		for(;i<t.end()-1;++i) {
			buf << *i << "<>";
			n += i->length() + 2;
		}
		buf << *i;
		n += i->length();
		return n;
	}
};

struct Blast_tab_format : public Output_format
{

	Blast_tab_format()
	{ }

	virtual void print_match(const DAA_match_record &r, Text_buffer &out) const
	{
		out << r.query_name() << '\t';

		if(r.subject_name.find_first_of('\1') == string::npos)
			out << r.subject_name;
		else
			this->print_salltitles(out, r.subject_name.c_str());

		out		<< '\t'
				<< (double)r.identities*100/r.len << '\t'
				<< r.len << '\t'
				<< r.mismatches << '\t'
				<< r.gap_openings << '\t'
				<< r.query_begin+1 << '\t'
				<< r.query_end()+1 << '\t'
				<< r.subject_begin+1 << '\t'
				<< r.subject_begin+r.subject_len << '\t';
		out.print_e(score_matrix.evalue(r.score, (size_t)r.db_letters(), (unsigned)r.query().size()));
		out << '\t' << score_matrix.bitscore(r.score) << '\n';
	}

	virtual ~Blast_tab_format()
	{ }

};

struct Sam_format : public Output_format
{

	Sam_format()
	{ }

	virtual void print_match(const DAA_match_record &r, Text_buffer &out) const
	{
		out << r.query_name() << '\t'
				<< '0' << '\t';

		if(r.subject_name.find_first_of('\1') == string::npos)
			out << r.subject_name;
		else
			this->print_salltitles(out, r.subject_name.c_str());

		out		<< '\t'
				<< r.subject_begin+1 << '\t'
				<< "255" << '\t';

		print_cigar(r, out);

		out << '\t'
				<< '*' << '\t'
				<< '0' << '\t'
				<< '0' << '\t'
				<< sequence (&r.query()[r.translated_query_begin], r.translated_query_len) << '\t'
				<< '*' << '\t'
				<< "AS:i:" << (uint32_t)score_matrix.bitscore(r.score) << '\t'
				<< "NM:i:" << r.len - r.identities << '\t'
				<< "ZL:i:" << r.total_subject_len << '\t'
				<< "ZR:i:" << r.score << '\t'
				<< "ZE:f:";
		out.print_e(score_matrix.evalue(r.score, (size_t)r.db_letters(), (unsigned)r.query().size()));
		out << '\t'
				<< "ZI:i:" << r.identities*100/r.len << '\t'
				<< "ZF:i:" << blast_frame(r.frame) << '\t'
				<< "ZS:i:" << r.query_begin+1 << '\t'
				<< "MD:Z:";

		print_md(r, out);
		out << '\n';
	}

	void print_md(const DAA_match_record &r, Text_buffer &buf) const
	{
		unsigned matches = 0, del = 0;
		for(Packed_transcript::Const_iterator i = r.transcript.begin(); i.good(); ++i) {
			switch(i->op) {
			case op_match:
				del = 0;
				matches += i->count;
				break;
			case op_insertion:
				break;
			case op_substitution:
				if(matches > 0) {
					buf << matches;
					matches = 0;
				} else if(del > 0) {
					buf << '0';
					del = 0;
				}
				buf << value_traits.alphabet[(long)i->letter];
				break;
			case op_deletion:
				if(matches > 0) {
					buf << matches;
					matches = 0;
				}
				if(del == 0)
					buf << '^';
				buf << value_traits.alphabet[(long)i->letter];
				++del;
			}
		}
		if(matches > 0)
			buf << matches;
	}

	void print_cigar(const DAA_match_record &r, Text_buffer &buf) const
	{
		static const unsigned map[] = { 0, 1, 2, 0 };
		static const char letter[] = { 'M', 'I', 'D' };
		unsigned n = 0, op = 0;
		for(Packed_transcript::Const_iterator i = r.transcript.begin(); i.good(); ++i) {
			if(map[i->op] == op)
				n += i->count;
			else {
				if(n > 0)
					buf << n << letter[op];
				n = i->count;
				op = map[i->op];
			}
		}
		if(n > 0)
			buf << n << letter[op];
	}

	virtual void print_header(Compressed_ostream &f, int mode) const
	{
		static const char* mode_str[] = { 0, 0, "BlastP", "BlastX", "BlastN" };
		string line = string("@HD\tVN:1.5\tSO:query\n\
@PG\tPN:DIAMOND\n\
@mm\t") + mode_str[mode] + "\n\
@CO\t" + mode_str[mode] + "-like alignments\n\
@CO\tReporting AS: bitScore, ZR: rawScore, ZE: expected, ZI: percent identity, ZL: reference length, ZF: frame, ZS: query start DNA coordinate\n";
		f.write(line.c_str(), line.length());
	}

	virtual ~Sam_format()
	{ }

};

const Output_format& get_output_format()
{
	static const Sam_format sam;
	static const Blast_tab_format tab;
	if(config.output_format == "tab")
		return tab;
	else if(config.output_format == "sam")
		return sam;
	else
		throw std::runtime_error("Invalid output format. Allowed values: tab,sam");
	return tab;
}

#endif /* OUTPUT_FORMAT_H_ */

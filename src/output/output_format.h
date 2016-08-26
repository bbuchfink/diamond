/****
Copyright (c) 2014-2016, University of Tuebingen, Benjamin Buchfink
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

#ifndef OUTPUT_FORMAT_H_
#define OUTPUT_FORMAT_H_

#include <exception>
#include "../basic/match.h"
#include "../align/match_func.h"
#include "../output/daa_file.h"
#include "../output/daa_record.h"
#include "../util/compressed_stream.h"
#include "../basic/score_matrix.h"

struct Output_format
{
	Output_format(unsigned code):
		code(code)
	{}
	virtual void print_query_intro(size_t query_num, const char *query_name, unsigned query_len, Text_buffer &out) const
	{}
	virtual void print_query_epilog(Text_buffer &out) const
	{}
	virtual void print_match(const Hsp_context& r, Text_buffer &out) const
	{}
	virtual void print_header(Output_stream &f, int mode, const char *matrix, int gap_open, int gap_extend, double evalue, const char *first_query_name, unsigned first_query_len) const
	{ }
	virtual void print_footer(Output_stream &f) const
	{ }
	virtual ~Output_format()
	{ }
	static void print_salltitles(Text_buffer &buf, const char *id, bool full_titles, bool all_titles)
	{
		if (!all_titles) {
			buf.write_until(id, full_titles ? "\1" : Const::id_delimiters);
			return;
		}
		if (strchr(id, '\1') == 0) {
			buf.write_until(id, full_titles ? "\1" : Const::id_delimiters);
			return;
		}
		//size_t n = 0;
		const vector<string> t (tokenize(id, "\1"));
		vector<string>::const_iterator i=t.begin();
		for(;i<t.end()-1;++i) {
			if (full_titles)
				buf << *i << "<>";
			else {
				buf.write_until(i->c_str(), Const::id_delimiters);
				buf << ";";
			}
			//n += i->length() + 2;
		}
		if(full_titles)
			buf << *i;
		else
			buf.write_until(i->c_str(), Const::id_delimiters);
		//n += i->length();
		//return n;
	}
	operator unsigned() const
	{
		return code;
	}
	unsigned code;
	enum { daa, blast_tab, blast_xml, sam };
};

extern auto_ptr<Output_format> output_format;

struct DAA_format : public Output_format
{
	DAA_format():
		Output_format(daa)
	{}
};

struct Blast_tab_format : public Output_format
{
	static const char* field_str[];
	Blast_tab_format();
	virtual void print_match(const Hsp_context& r, Text_buffer &out) const;
	virtual ~Blast_tab_format()
	{ }
	vector<unsigned> fields;
};

struct Sam_format : public Output_format
{

	Sam_format():
		Output_format(sam)
	{ }

	virtual void print_match(const Hsp_context& r, Text_buffer &out) const
	{
		out.write_until(r.query_name, Const::id_delimiters);
		out << '\t' << '0' << '\t';

		const bool lt = (config.salltitles || (config.command == Config::view)) ? true : false;
		this->print_salltitles(out, r.subject_name, lt, lt);

		out << '\t'
			<< r.subject_range().begin_ + 1 << '\t'
			<< "255" << '\t';

		print_cigar(r, out);

		out << '\t'
			<< '*' << '\t'
			<< '0' << '\t'
			<< '0' << '\t'
			<< sequence(&r.query[r.query_range().begin_], r.query_range().length()) << '\t'
			<< '*' << '\t'
			<< "AS:i:" << (uint32_t)score_matrix.bitscore(r.score()) << '\t'
			<< "NM:i:" << r.length() - r.identities() << '\t'
			<< "ZL:i:" << r.subject_len << '\t'
			<< "ZR:i:" << r.score() << '\t'
			<< "ZE:f:";
		out.print_e(score_matrix.evalue(r.score(), config.db_size, (unsigned)r.query.length()));
		out << '\t'
			<< "ZI:i:" << r.identities() * 100 / r.length() << '\t'
			<< "ZF:i:" << blast_frame(r.frame()) << '\t'
			<< "ZS:i:" << r.oriented_query_range().begin_ + 1 << '\t'
			<< "MD:Z:";

		print_md(r, out);
		out << '\n';
	}

	void print_md(const Hsp_context &r, Text_buffer &buf) const
	{
		unsigned matches = 0, del = 0;
		for(Packed_transcript::Const_iterator i = r.begin_old(); i.good(); ++i) {
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

	void print_cigar(const Hsp_context &r, Text_buffer &buf) const
	{
		static const unsigned map[] = { 0, 1, 2, 0 };
		static const char letter[] = { 'M', 'I', 'D' };
		unsigned n = 0, op = 0;
		for(Packed_transcript::Const_iterator i = r.begin_old(); i.good(); ++i) {
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

	virtual void print_header(Output_stream &f, int mode, const char *matrix, int gap_open, int gap_extend, double evalue, const char *first_query_name, unsigned first_query_len) const
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

struct XML_format : public Output_format
{
	XML_format():
		Output_format(blast_xml)
	{}
	virtual void print_match(const Hsp_context &r, Text_buffer &out) const;
	virtual void print_header(Output_stream &f, int mode, const char *matrix, int gap_open, int gap_extend, double evalue, const char *first_query_name, unsigned first_query_len) const;
	virtual void print_query_intro(size_t query_num, const char *query_name, unsigned query_len, Text_buffer &out) const;
	virtual void print_query_epilog(Text_buffer &out) const;
	virtual void print_footer(Output_stream &f) const;
	virtual ~XML_format()
	{ }
};

Output_format* get_output_format();

#endif /* OUTPUT_FORMAT_H_ */

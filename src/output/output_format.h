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
	virtual void print_query_intro(size_t query_num, const char *query_name, unsigned query_len, Text_buffer &out, bool unaligned) const
	{}
	virtual void print_query_epilog(Text_buffer &out, bool unaligned) const
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
	enum { daa, blast_tab, blast_xml, sam, blast_pairwise };
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
	virtual void print_query_intro(size_t query_num, const char *query_name, unsigned query_len, Text_buffer &out, bool unaligned) const;
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
	virtual void print_match(const Hsp_context& r, Text_buffer &out) const;
	virtual void print_header(Output_stream &f, int mode, const char *matrix, int gap_open, int gap_extend, double evalue, const char *first_query_name, unsigned first_query_len) const;
	virtual void print_query_intro(size_t query_num, const char *query_name, unsigned query_len, Text_buffer &out, bool unaligned) const;
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
	virtual void print_query_intro(size_t query_num, const char *query_name, unsigned query_len, Text_buffer &out, bool unaligned) const;
	virtual void print_query_epilog(Text_buffer &out, bool unaligned) const;
	virtual void print_footer(Output_stream &f) const;
	virtual ~XML_format()
	{ }
};

struct Pairwise_format : public Output_format
{
	Pairwise_format() :
		Output_format(blast_pairwise)
	{}
	virtual void print_match(const Hsp_context &r, Text_buffer &out) const;
	virtual void print_header(Output_stream &f, int mode, const char *matrix, int gap_open, int gap_extend, double evalue, const char *first_query_name, unsigned first_query_len) const;
	virtual void print_query_intro(size_t query_num, const char *query_name, unsigned query_len, Text_buffer &out, bool unaligned) const;
	virtual void print_query_epilog(Text_buffer &out, bool unaligned) const;
	virtual void print_footer(Output_stream &f) const;
	virtual ~Pairwise_format()
	{ }
};

Output_format* get_output_format();

#endif /* OUTPUT_FORMAT_H_ */

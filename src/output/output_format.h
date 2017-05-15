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
	static void print_title(Text_buffer &buf, const char *id, bool full_titles, bool all_titles, const char *separator);
	operator unsigned() const
	{
		return code;
	}
	unsigned code;
	enum { daa, blast_tab, blast_xml, sam, blast_pairwise, null };
};

extern auto_ptr<Output_format> output_format;

struct Null_format : public Output_format
{
	Null_format() :
		Output_format(null)
	{}
};

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
	{
		config.salltitles = true;
	}
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
	{
		config.salltitles = true;
	}
	virtual void print_match(const Hsp_context &r, Text_buffer &out) const;
	virtual void print_header(Output_stream &f, int mode, const char *matrix, int gap_open, int gap_extend, double evalue, const char *first_query_name, unsigned first_query_len) const;
	virtual void print_query_intro(size_t query_num, const char *query_name, unsigned query_len, Text_buffer &out, bool unaligned) const;
	virtual void print_query_epilog(Text_buffer &out, bool unaligned) const;
	virtual void print_footer(Output_stream &f) const;
	virtual ~Pairwise_format()
	{ }
};

Output_format* get_output_format();
void print_hsp(Hsp_data &hsp, sequence query);

#endif /* OUTPUT_FORMAT_H_ */

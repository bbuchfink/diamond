/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen
						
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

#pragma once
#include <exception>
#include <memory>
#include "../basic/match.h"
#include "../output/daa_file.h"
#include "../output/daa_record.h"
#include "../util/io/output_file.h"
#include "../basic/score_matrix.h"
#include "../util/escape_sequences.h"
#include "../basic/parameters.h"
#include "../data/metadata.h"
#include "../util/io/consumer.h"

struct Output_format
{
	Output_format(unsigned code):
		code(code),
		needs_taxon_id_lists(false),
		needs_taxon_nodes(false),
		needs_taxon_scientific_names(false),
		needs_taxon_ranks(false),
		needs_transcript(true),
		needs_stats(false)
	{}
	virtual void print_query_intro(size_t query_num, const char *query_name, unsigned query_len, TextBuffer &out, bool unaligned) const
	{}
	virtual void print_query_epilog(TextBuffer &out, const char *query_title, bool unaligned, const Parameters &parameters) const
	{}
	virtual void print_match(const Hsp_context& r, const Metadata &metadata, TextBuffer &out)
	{}
	virtual void print_header(Consumer &f, int mode, const char *matrix, int gap_open, int gap_extend, double evalue, const char *first_query_name, unsigned first_query_len) const
	{ }
	virtual void print_footer(Consumer &f) const
	{ }
	virtual Output_format* clone() const = 0;
	virtual ~Output_format()
	{ }
	static void print_title(TextBuffer &buf, const char *id, bool full_titles, bool all_titles, const char *separator, const EscapeSequences *esc = 0);
	operator unsigned() const
	{
		return code;
	}
	unsigned code;
	bool needs_taxon_id_lists, needs_taxon_nodes, needs_taxon_scientific_names, needs_taxon_ranks, needs_transcript, needs_stats;
	enum { daa, blast_tab, blast_xml, sam, blast_pairwise, null, taxon, paf, bin1 };
};

extern std::unique_ptr<Output_format> output_format;

struct Null_format : public Output_format
{
	Null_format() :
		Output_format(null)
	{}
	virtual Output_format* clone() const
	{
		return new Null_format(*this);
	}
};

struct DAA_format : public Output_format
{
	DAA_format():
		Output_format(daa)
	{}
	virtual Output_format* clone() const
	{
		return new DAA_format(*this);
	}
};

struct OutputField {
	const std::string key, description;
	const bool need_transcript, need_stats;
};

struct Blast_tab_format : public Output_format
{
	static const std::vector<OutputField> field_def;
	Blast_tab_format();
	virtual void print_header(Consumer &f, int mode, const char *matrix, int gap_open, int gap_extend, double evalue, const char *first_query_name, unsigned first_query_len) const override;
	virtual void print_query_intro(size_t query_num, const char *query_name, unsigned query_len, TextBuffer &out, bool unaligned) const override;
	virtual void print_match(const Hsp_context& r, const Metadata &metadata, TextBuffer &out) override;
	virtual ~Blast_tab_format()
	{ }
	virtual Output_format* clone() const override
	{
		return new Blast_tab_format(*this);
	}
	vector<unsigned> fields;
};

struct PAF_format : public Output_format
{
	PAF_format():
		Output_format(paf)
	{}
	virtual void print_query_intro(size_t query_num, const char *query_name, unsigned query_len, TextBuffer &out, bool unaligned) const;
	virtual void print_match(const Hsp_context& r, const Metadata &metadata, TextBuffer &out);
	virtual ~PAF_format()
	{ }
	virtual Output_format* clone() const
	{
		return new PAF_format(*this);
	}
};

struct Sam_format : public Output_format
{
	Sam_format():
		Output_format(sam)
	{ }
	virtual void print_match(const Hsp_context& r, const Metadata &metadata, TextBuffer &out) override;
	virtual void print_header(Consumer &f, int mode, const char *matrix, int gap_open, int gap_extend, double evalue, const char *first_query_name, unsigned first_query_len) const override;
	virtual void print_query_intro(size_t query_num, const char *query_name, unsigned query_len, TextBuffer &out, bool unaligned) const override;
	virtual ~Sam_format()
	{ }
	virtual Output_format* clone() const override
	{
		return new Sam_format(*this);
	}
};

struct XML_format : public Output_format
{
	XML_format():
		Output_format(blast_xml)
	{
		config.salltitles = true;
	}
	virtual void print_match(const Hsp_context& r, const Metadata &metadata, TextBuffer &out) override;
	virtual void print_header(Consumer &f, int mode, const char *matrix, int gap_open, int gap_extend, double evalue, const char *first_query_name, unsigned first_query_len) const override;
	virtual void print_query_intro(size_t query_num, const char *query_name, unsigned query_len, TextBuffer &out, bool unaligned) const override;
	virtual void print_query_epilog(TextBuffer &out, const char *query_title, bool unaligned, const Parameters &parameters) const override;
	virtual void print_footer(Consumer &f) const override;
	virtual ~XML_format()
	{ }
	virtual Output_format* clone() const override
	{
		return new XML_format(*this);
	}
};

struct Pairwise_format : public Output_format
{
	Pairwise_format() :
		Output_format(blast_pairwise)
	{
		config.salltitles = true;
	}
	virtual void print_match(const Hsp_context& r, const Metadata &metadata, TextBuffer &out) override;
	virtual void print_header(Consumer &f, int mode, const char *matrix, int gap_open, int gap_extend, double evalue, const char *first_query_name, unsigned first_query_len) const override;
	virtual void print_query_intro(size_t query_num, const char *query_name, unsigned query_len, TextBuffer &out, bool unaligned) const override;
	virtual void print_query_epilog(TextBuffer &out, const char *query_title, bool unaligned, const Parameters &parameters) const override;
	virtual void print_footer(Consumer &f) const override;
	virtual ~Pairwise_format()
	{ }
	virtual Output_format* clone() const override
	{
		return new Pairwise_format(*this);
	}
};

struct Taxon_format : public Output_format
{
	Taxon_format() :
		Output_format(taxon),
		taxid(0),
		evalue(std::numeric_limits<double>::max())
	{
		needs_taxon_id_lists = true;
		needs_taxon_nodes = true;
	}
	virtual void print_match(const Hsp_context& r, const Metadata &metadata, TextBuffer &out) override;
	virtual void print_query_epilog(TextBuffer &out, const char *query_title, bool unaligned, const Parameters &parameters) const override;
	virtual ~Taxon_format()
	{ }
	virtual Output_format* clone() const override
	{
		return new Taxon_format(*this);
	}
	unsigned taxid;
	double evalue;
};

struct Bin1_format : public Output_format
{
	Bin1_format():
		Output_format(bin1)
	{}
	virtual void print_query_intro(size_t query_num, const char *query_name, unsigned query_len, TextBuffer &out, bool unaligned) const override;
	virtual void print_match(const Hsp_context& r, const Metadata &metadata, TextBuffer &out) override;
	virtual ~Bin1_format()
	{ }
	virtual Output_format* clone() const override
	{
		return new Bin1_format(*this);
	}
};

struct Clustering_format : public Output_format
{
	const string* const format;
	Clustering_format(const string* const format): Output_format(bin1), format(format) {
	}
	virtual void print_query_intro(size_t query_num, const char *query_name, unsigned query_len, TextBuffer &out, bool unaligned) const override;
	virtual void print_match(const Hsp_context& r, const Metadata &metadata, TextBuffer &out) override;
	virtual ~Clustering_format()
	{ }
	virtual Output_format* clone() const override
	{
		return new Clustering_format(*this);
	}
};

Output_format* get_output_format();
void init_output(bool have_taxon_id_lists, bool have_taxon_nodes, bool have_taxon_scientific_names);
void print_hsp(Hsp &hsp, const TranslatedSequence &query);
void print_cigar(const Hsp_context &r, TextBuffer &buf);


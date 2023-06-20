/****
DIAMOND protein aligner
Copyright (C) 2013-2021 Max Planck Society for the Advancement of Science e.V.
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
#include <algorithm>
#include <exception>
#include <memory>
#include "../basic/match.h"
#include "../output/daa/daa_file.h"
#include "../output/daa/daa_record.h"
#include "../util/io/output_file.h"
#include "../stats/score_matrix.h"
#include "../util/escape_sequences.h"
#include "../util/io/consumer.h"
#include "../util/enum.h"
#include "../run/config.h"
#include "../dp/flags.h"
#include "../data/taxonomy.h"
#include "def.h"

namespace Output {

struct Info {	
	SeqInfo query;
	bool unaligned;
	SequenceFile* db;
	TextBuffer& out;
	Extension::Stats stats;
	AccessionParsing acc_stats;
};

}

struct OutputFormat
{
	OutputFormat(unsigned code, HspValues hsp_values = HspValues::TRANSCRIPT, Output::Flags flags = Output::Flags::NONE):
		code(code),
		needs_taxon_id_lists(false),
		needs_taxon_nodes(false),
		needs_taxon_scientific_names(false),
		needs_taxon_ranks(false),
		needs_paired_end_info(false),
		hsp_values(hsp_values),
		flags(flags)
	{}
	virtual void print_query_intro(Output::Info& info) const
	{}
	virtual void print_query_epilog(Output::Info& info) const
	{}
	virtual void print_match(const HspContext& r, Output::Info& info)
	{}
	virtual void print_header(Consumer &f, int mode, const char *matrix, int gap_open, int gap_extend, double evalue, const char *first_query_name, unsigned first_query_len) const
	{ }
	virtual void print_footer(Consumer &f) const
	{ }
	virtual OutputFormat* clone() const = 0;
	virtual ~OutputFormat()
	{ }
	static void print_title(TextBuffer &buf, const char *id, bool full_titles, bool all_titles, const char *separator, const EscapeSequences *esc = 0,  bool json_array = false);
	operator unsigned() const
	{
		return code;
	}
	unsigned code;
	bool needs_taxon_id_lists, needs_taxon_nodes, needs_taxon_scientific_names, needs_taxon_ranks, needs_paired_end_info;
	HspValues hsp_values;
	Output::Flags flags;
	enum { daa, blast_tab, blast_xml, sam, blast_pairwise, null, taxon, paf, bin1, EDGE ,json};
};

struct Null_format : public OutputFormat
{
	Null_format() :
		OutputFormat(null)
	{}
	virtual OutputFormat* clone() const
	{
		return new Null_format(*this);
	}
};

struct DAA_format : public OutputFormat
{
	DAA_format();
	virtual OutputFormat* clone() const
	{
		return new DAA_format(*this);
	}
};

struct OutputField {
	const std::string key, clust_key, description;
	const HspValues hsp_values;
	const Output::Flags flags;
};

enum class Header { NONE, SIMPLE, VERBOSE };

struct Blast_tab_format : public OutputFormat
{
    Blast_tab_format(bool json = false);
    static const std::vector<OutputField> field_def;
    virtual void print_header(Consumer &f, int mode, const char *matrix, int gap_open, int gap_extend, double evalue, const char *first_query_name, unsigned first_query_len) const override;
    virtual void print_footer(Consumer &f) const override;
    virtual void print_query_intro(Output::Info& info) const override;
    virtual void print_match(const HspContext& r, Output::Info& info) override;
    virtual ~Blast_tab_format()
    { }
    virtual OutputFormat* clone() const override
    {
        return new Blast_tab_format(*this);
    }
    static Header header_format(unsigned workflow);
    void output_header(Consumer& f, bool cluster) const;
    std::vector<int64_t> fields;
    bool is_json;
};


struct PAF_format : public OutputFormat
{
	PAF_format():
		OutputFormat(paf, HspValues::TRANSCRIPT, Output::Flags::NONE)
	{}
	virtual void print_query_intro(Output::Info& info) const override;
	virtual void print_match(const HspContext& r, Output::Info& info) override;
	virtual ~PAF_format()
	{ }
	virtual OutputFormat* clone() const
	{
		return new PAF_format(*this);
	}
};

struct Sam_format : public OutputFormat
{
	Sam_format():
		OutputFormat(sam, HspValues::TRANSCRIPT, Output::Flags::NONE)
	{ }
	virtual void print_match(const HspContext& r, Output::Info& info) override;
	virtual void print_header(Consumer &f, int mode, const char *matrix, int gap_open, int gap_extend, double evalue, const char *first_query_name, unsigned first_query_len) const override;
	virtual void print_query_intro(Output::Info& info) const override;
	virtual ~Sam_format()
	{ }
	virtual OutputFormat* clone() const override
	{
		return new Sam_format(*this);
	}
};

struct XML_format : public OutputFormat
{
	XML_format():
		OutputFormat(blast_xml, HspValues::TRANSCRIPT, Output::Flags::FULL_TITLES)
	{
		config.salltitles = true;
	} 
	virtual void print_match(const HspContext& r, Output::Info& info) override;
	virtual void print_header(Consumer &f, int mode, const char *matrix, int gap_open, int gap_extend, double evalue, const char *first_query_name, unsigned first_query_len) const override;
	virtual void print_query_intro(Output::Info& info) const override;
	virtual void print_query_epilog(Output::Info& info) const override;
	virtual void print_footer(Consumer &f) const override;
	virtual ~XML_format()
	{ }
	virtual OutputFormat* clone() const override
	{
		return new XML_format(*this);
	}
};

struct Pairwise_format : public OutputFormat
{
	Pairwise_format() :
		OutputFormat(blast_pairwise, HspValues::TRANSCRIPT, Output::Flags::FULL_TITLES)
	{
		config.salltitles = true;
	}
	virtual void print_match(const HspContext& r, Output::Info& info) override;
	virtual void print_header(Consumer &f, int mode, const char *matrix, int gap_open, int gap_extend, double evalue, const char *first_query_name, unsigned first_query_len) const override;
	virtual void print_query_intro(Output::Info& info) const override;
	virtual void print_query_epilog(Output::Info& infos) const override;
	virtual void print_footer(Consumer &f) const override;
	virtual ~Pairwise_format()
	{ }
	virtual OutputFormat* clone() const override
	{
		return new Pairwise_format(*this);
	}
};

struct Taxon_format : public OutputFormat
{
	Taxon_format() :
		OutputFormat(taxon, HspValues::NONE, Output::Flags::NONE),
		taxid(0),
		evalue(std::numeric_limits<double>::max())
	{
		needs_taxon_id_lists = true;
		needs_taxon_nodes = true;
	}
	virtual void print_match(const HspContext& r, Output::Info& info) override;
	virtual void print_query_epilog(Output::Info& info) const override;
	virtual ~Taxon_format()
	{ }
	virtual OutputFormat* clone() const override
	{
		return new Taxon_format(*this);
	}
	unsigned taxid;
	double evalue;
};

struct Bin1_format : public OutputFormat
{
	Bin1_format():
		OutputFormat(bin1)
	{}
	virtual void print_query_intro(Output::Info& info) const override;
	virtual void print_match(const HspContext& r, Output::Info& info) override;
	virtual ~Bin1_format()
	{ }
	virtual OutputFormat* clone() const override
	{
		return new Bin1_format(*this);
	}
};

struct Clustering_format : public OutputFormat
{
	std::string format;
	Clustering_format(const std::string* const format);
	virtual void print_match(const HspContext& r, Output::Info& info) override;
	virtual ~Clustering_format()
	{ }
	virtual OutputFormat* clone() const override
	{
		return new Clustering_format(*this);
	}
};

namespace Output { namespace Format {

struct Edge : public OutputFormat
{
	struct Data {
		OId query, target;
		float qcovhsp, scovhsp;
		double evalue;
	};
	Edge() :
		OutputFormat(EDGE, HspValues::COORDS)
	{}
	virtual void print_match(const HspContext& r, Output::Info& info) override;
	virtual ~Edge()
	{ }
	virtual OutputFormat* clone() const override
	{
		return new Edge(*this);
	}
};

}}

OutputFormat* get_output_format();
OutputFormat* init_output(const int64_t max_target_seqs);
void print_hsp(Hsp &hsp, const TranslatedSequence &query);
void print_cigar(const HspContext &r, TextBuffer &buf);

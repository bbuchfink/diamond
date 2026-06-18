/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once
#include <functional>
#include <map>
#include "basic/match.h"
#include "dp/flags.h"
#include "def.h"
#include "util/sequence/sequence.h"
#include "util/escape_sequences.h"
#include "basic/config.h"
#include "search/seed_array/flags.h"
#include "util/io/file.h"
#include "util/optional.h"

struct SequenceFile;

namespace Output {

struct Info {
	SeqInfo query;
	bool unaligned;
	SequenceFile* db;
	TextBuffer& out;
	Util::Seq::AccessionParsing acc_stats;
	optional<uint64_t> db_seqs, db_letters;
};

}

enum class FieldId : unsigned {
	QSeqId = 0,
	QLen = 4,
	SSeqId = 5,
	SAllSeqId = 6,
	SLen = 12,
	QStart = 13,
	QEnd = 14,
	SStart = 15,
	SEnd = 16,
	QSeq = 17,
	SSeq = 18,
	EValue = 19,
	BitScore = 20,
	Score = 21,
	Length = 22,
	PIdent = 23,
	NIdent = 24,
	Mismatch = 25,
	Positive = 26,
	GapOpen = 27,
	Gaps = 28,
	PPos = 29,
	QFrame = 31,
	BTop = 33,
	STaxIds = 34,
	SSciNames = 35,
	SSKingdoms = 38,
	STitle = 39,
	SAllTitles = 40,
	SStrand = 41,
	QCovHsp = 43,
	QTitle = 45,
	FullSSeq = 48,
	QQual = 49,
	QNum = 50,
	SNum = 51,
	SCovHsp = 52,
	FullQQual = 53,
	FullQSeq = 54,
	QSeqGapped = 55,
	SSeqGapped = 56,
	QStrand = 57,
	Cigar = 58,
	SKingdoms = 59,
	SPhylums = 60,
	FullQSeqMate = 62,
	QSeqTranslated = 63,
	HspNum = 66,
	NormalizedBitscore = 67,
	NORMALIZED_NIDENT = 68,
	ApproxPIdent = 71,
	CorrectedBitScore = 72,
	NegEValue = 73,
	Reserved1 = 74,
	Reserved2 = 75,
	SLineages = 76,
	NormalizedBitscoreQuery = 77,
	COUNT = 78
};

struct OutputFormat
{
	OutputFormat(unsigned code, HspValues hsp_values = HspValues::TRANSCRIPT, Output::Flags flags = Output::Flags::NONE, char query_separator = '\0') :
		code(code),
		query_separator(query_separator),
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
	virtual void print_header(File &f, int mode, const char *matrix, int gap_open, int gap_extend, double evalue, const char *first_query_name, unsigned first_query_len) const
	{ }
	virtual void print_footer(File &f) const
	{ }
	virtual OutputFormat* clone() const = 0;
	virtual ~OutputFormat()
	{ }
	static void print_title(TextBuffer &buf, const char *id, bool full_titles, bool all_titles, const char *separator, const EscapeSequences *esc = 0,  bool json_array = false);
	operator unsigned() const
	{
		return code;
	}
	bool report_unaligned() const {
		return config.report_unaligned != 0 && (config.report_unaligned != -1 || flag_any(flags, Output::Flags::DEFAULT_REPORT_UNALIGNED));
	}
	const unsigned code;
	const char query_separator;
	bool needs_taxon_id_lists, needs_taxon_nodes, needs_taxon_scientific_names, needs_taxon_ranks, needs_paired_end_info;
	HspValues hsp_values;
	Output::Flags flags;
	enum { daa, blast_tab, blast_xml, sam, blast_pairwise, null, taxon, paf, bin1, EDGE, json };
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

struct DAAFormat : public OutputFormat
{
	DAAFormat();
	virtual OutputFormat* clone() const
	{
		return new DAAFormat(*this);
	}
};

enum class FieldId : unsigned;

struct OutputField {
	FieldId id;
	std::string key, clust_key, description;
	HspValues hsp_values;
	Output::Flags flags;
};

enum class Header { NONE, SIMPLE, VERBOSE };

struct TabularFormat;

struct FieldCallbacks {
	std::function<void(const TabularFormat&, const HspContext&, Output::Info&)> match;
	std::function<void(const TabularFormat&, Output::Info&)> query_intro;
};

struct TabularFormat : public OutputFormat
{
	TabularFormat(bool json = false);
	static std::map<FieldId, OutputField> field_def;
	static std::map<FieldId, FieldCallbacks> field_callbacks;
	virtual void print_header(File& f, int mode, const char* matrix, int gap_open, int gap_extend, double evalue, const char* first_query_name, unsigned first_query_len) const override;
	virtual void print_footer(File& f) const override;
	virtual void print_query_intro(Output::Info& info) const override;
	virtual void print_match(const HspContext& r, Output::Info& info) override;
	virtual ~TabularFormat()
	{
	}
	virtual OutputFormat* clone() const override
	{
		return new TabularFormat(*this);
	}
	static Header header_format(unsigned workflow);
	void output_header(File& f, bool cluster) const;
	std::vector<FieldId> fields;
	bool is_json;
};


struct PAFFormat : public OutputFormat
{
	PAFFormat():
		OutputFormat(paf, HspValues::TRANSCRIPT, Output::Flags::SSEQID | Output::Flags::DEFAULT_REPORT_UNALIGNED)
	{}
	virtual void print_query_intro(Output::Info& info) const override;
	virtual void print_match(const HspContext& r, Output::Info& info) override;
	virtual ~PAFFormat()
	{ }
	virtual OutputFormat* clone() const override
	{
		return new PAFFormat(*this);
	}
};

struct SamFormat : public OutputFormat
{
	SamFormat():
		OutputFormat(sam, HspValues::TRANSCRIPT, Output::Flags::SSEQID | Output::Flags::DEFAULT_REPORT_UNALIGNED)
	{ }
	virtual void print_match(const HspContext& r, Output::Info& info) override;
	virtual void print_header(File&f, int mode, const char *matrix, int gap_open, int gap_extend, double evalue, const char *first_query_name, unsigned first_query_len) const override;
	virtual void print_query_intro(Output::Info& info) const override;
	virtual ~SamFormat()
	{ }
	virtual OutputFormat* clone() const override
	{
		return new SamFormat(*this);
	}
};

struct XMLFormat : public OutputFormat
{
	XMLFormat():
		OutputFormat(blast_xml, HspValues::TRANSCRIPT, Output::Flags::FULL_TITLES | Output::Flags::SSEQID | Output::Flags::DEFAULT_REPORT_UNALIGNED | Output::Flags::ALL_SEQIDS)
	{} 
	virtual void print_match(const HspContext& r, Output::Info& info) override;
	virtual void print_header(File &f, int mode, const char *matrix, int gap_open, int gap_extend, double evalue, const char *first_query_name, unsigned first_query_len) const override;
	virtual void print_query_intro(Output::Info& info) const override;
	virtual void print_query_epilog(Output::Info& info) const override;
	virtual void print_footer(File&f) const override;
	virtual ~XMLFormat()
	{ }
	virtual OutputFormat* clone() const override
	{
		return new XMLFormat(*this);
	}
};

struct PairwiseFormat : public OutputFormat
{
	PairwiseFormat() :
		OutputFormat(blast_pairwise, HspValues::TRANSCRIPT, Output::Flags::FULL_TITLES | Output::Flags::SSEQID | Output::Flags::DEFAULT_REPORT_UNALIGNED | Output::Flags::ALL_SEQIDS)
	{}
	virtual void print_match(const HspContext& r, Output::Info& info) override;
	virtual void print_header(File&f, int mode, const char *matrix, int gap_open, int gap_extend, double evalue, const char *first_query_name, unsigned first_query_len) const override;
	virtual void print_query_intro(Output::Info& info) const override;
	virtual void print_query_epilog(Output::Info& infos) const override;
	virtual void print_footer(File&f) const override;
	virtual ~PairwiseFormat()
	{ }
	virtual OutputFormat* clone() const override
	{
		return new PairwiseFormat(*this);
	}
};

struct TaxonFormat : public OutputFormat
{
	TaxonFormat();
	virtual void print_match(const HspContext& r, Output::Info& info) override;
	virtual void print_query_epilog(Output::Info& info) const override;
	virtual ~TaxonFormat()
	{ }
	virtual OutputFormat* clone() const override
	{
		return new TaxonFormat(*this);
	}
	TaxId taxid;
	double evalue;
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
		OutputFormat(EDGE, HspValues::COORDS, Output::Flags::SSEQID)
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
OutputFormat* init_output(int64_t& max_target_seqs);
void print_hsp(Hsp &hsp, const TranslatedSequence &query);
void print_cigar(const HspContext &r, TextBuffer &buf);

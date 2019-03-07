/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

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
#include <set>
#include <numeric>
#include "../basic/match.h"
#include "output_format.h"
#include "../data/taxonomy.h"
#include "../data/queries.h"

using namespace std;

const char* Blast_tab_format::field_str[] = {
	"qseqid",		// 0 means Query Seq - id
	"qgi",			// 1 means Query GI
	"qacc",			// 2 means Query accesion
	"qaccver",		// 3 means Query accesion.version
	"qlen",			// 4 means Query sequence length
	"sseqid",		// 5 means Subject Seq - id
	"sallseqid",	// 6 means All subject Seq - id(s), separated by a ';'
	"sgi",			// 7 means Subject GI
	"sallgi",		// 8 means All subject GIs
	"sacc",			// 9 means Subject accession
	"saccver",		// 10 means Subject accession.version
	"sallacc",		// 11 means All subject accessions
	"slen",			// 12 means Subject sequence length
	"qstart",		// 13 means Start of alignment in query
	"qend",			// 14 means End of alignment in query
	"sstart",		// 15 means Start of alignment in subject
	"send",			// 16 means End of alignment in subject
	"qseq",			// 17 means Aligned part of query sequence
	"sseq",			// 18 means Aligned part of subject sequence
	"evalue",		// 19 means Expect value
	"bitscore",		// 20 means Bit score
	"score",		// 21 means Raw score
	"length",		// 22 means Alignment length
	"pident",		// 23 means Percentage of identical matches
	"nident",		// 24 means Number of identical matches
	"mismatch",		// 25 means Number of mismatches
	"positive",		// 26 means Number of positive - scoring matches
	"gapopen",		// 27 means Number of gap openings
	"gaps",			// 28 means Total number of gaps
	"ppos",			// 29 means Percentage of positive - scoring matches
	"frames",		// 30 means Query and subject frames separated by a '/'
	"qframe",		// 31 means Query frame
	"sframe",		// 32 means Subject frame
	"btop",			// 33 means Blast traceback operations(BTOP)
	"staxids",		// 34 means unique Subject Taxonomy ID(s), separated by a ';'	(in numerical order)
	"sscinames",	// 35 means unique Subject Scientific Name(s), separated by a ';'
	"scomnames",	// 36 means unique Subject Common Name(s), separated by a ';'
	"sblastnames",	// 37 means unique Subject Blast Name(s), separated by a ';'	(in alphabetical order)
	"sskingdoms",	// 38 means unique Subject Super Kingdom(s), separated by a ';'	(in alphabetical order)
	"stitle",		// 39 means Subject Title
	"salltitles",	// 40 means All Subject Title(s), separated by a '<>'
	"sstrand",		// 41 means Subject Strand
	"qcovs",		// 42 means Query Coverage Per Subject
	"qcovhsp",		// 43 means Query Coverage Per HSP
	"qcovus",		// 44 means Query Coverage Per Unique Subject(blastn only)
	"qtitle",		// 45 means Query title
	"swdiff",		// 46
	"time", 		// 47
	"full_sseq",	// 48
	"qqual",		// 49
	"qnum",			// 50
	"snum",			// 51
	"scovhsp",		// 52
	"full_qqual",	// 53
	"full_qseq"		// 54
};

const char* Blast_tab_format::field_desc[] = {
	"Query ID",		// 0 means Query Seq - id
	"qgi",			// 1 means Query GI
	"qacc",			// 2 means Query accesion
	"qaccver",		// 3 means Query accesion.version
	"Query length",	// 4 means Query sequence length
	"Subject ID",		// 5 means Subject Seq - id
	"Subject IDs",	// 6 means All subject Seq - id(s), separated by a ';'
	"sgi",			// 7 means Subject GI
	"sallgi",		// 8 means All subject GIs
	"sacc",			// 9 means Subject accession
	"saccver",		// 10 means Subject accession.version
	"sallacc",		// 11 means All subject accessions
	"Subject length",			// 12 means Subject sequence length
	"Start of alignment in query",		// 13 means Start of alignment in query
	"End of alignment in query",			// 14 means End of alignment in query
	"Start of alignment in subject",		// 15 means Start of alignment in subject
	"End of alignment in subject",			// 16 means End of alignment in subject
	"Aligned part of query sequence",			// 17 means Aligned part of query sequence
	"Aligned part of subject sequence",			// 18 means Aligned part of subject sequence
	"Expected value",		// 19 means Expect value
	"Bit score",		// 20 means Bit score
	"Raw score",		// 21 means Raw score
	"Alignment length",		// 22 means Alignment length
	"Percentage of identical matches",		// 23 means Percentage of identical matches
	"Number of identical matches",		// 24 means Number of identical matches
	"Number of mismatches",		// 25 means Number of mismatches
	"Number of positive-scoring matches",		// 26 means Number of positive - scoring matches
	"Number of gap openings",		// 27 means Number of gap openings
	"Total number of gaps",			// 28 means Total number of gaps
	"Percentage of positive-scoring matches",			// 29 means Percentage of positive - scoring matches
	"frames",		// 30 means Query and subject frames separated by a '/'
	"Query frame",		// 31 means Query frame
	"sframe",		// 32 means Subject frame
	"Blast traceback operations",			// 33 means Blast traceback operations(BTOP)
	"Subject Taxonomy IDs",		// 34 means unique Subject Taxonomy ID(s), separated by a ';'	(in numerical order)
	"Subject scientific names",	// 35 means unique Subject Scientific Name(s), separated by a ';'
	"scomnames",	// 36 means unique Subject Common Name(s), separated by a ';'
	"sblastnames",	// 37 means unique Subject Blast Name(s), separated by a ';'	(in alphabetical order)
	"sskingdoms",	// 38 means unique Subject Super Kingdom(s), separated by a ';'	(in alphabetical order)
	"Subject title",		// 39 means Subject Title
	"Subject titles",	// 40 means All Subject Title(s), separated by a '<>'
	"sstrand",		// 41 means Subject Strand
	"qcovs",		// 42 means Query Coverage Per Subject
	"Query coverage per HSP",		// 43 means Query Coverage Per HSP
	"qcovus",		// 44 means Query Coverage Per Unique Subject(blastn only)
	"Query title",		// 45 means Query title
	"swdiff",		// 46
	"time", 		// 47
	"Subject sequence",	// 48
	"Aligned part of query quality values",		// 49
	"qnum",			// 50
	"snum",			// 51
	"scovhsp",		// 52
	"Query quality values",	// 53
	"Query sequence"		// 54
};

Blast_tab_format::Blast_tab_format() :
	Output_format(blast_tab)
{
	static const unsigned stdf[] = { 0, 5, 23, 22, 25, 27, 13, 14, 15, 16, 19, 20 };
	const vector<string> &f = config.output_format;
	if (f.size() <= 1) {
		fields = vector<unsigned>(stdf, stdf + 12);
		return;
	}
	for (vector<string>::const_iterator i = f.begin() + 1; i != f.end(); ++i) {
		int j = get_idx(field_str, sizeof(field_str) / sizeof(field_str[0]), i->c_str());
		if(j == -1)
			throw std::runtime_error(string("Invalid output field: ") + *i);
		if (j == 34)
			needs_taxon_id_lists = true;
		if (j == 35) {
			needs_taxon_scientific_names = true;
			needs_taxon_id_lists = true;
		}
		fields.push_back(j);
		if (j == 6 || j == 39 || j == 40 || j == 34)
			config.salltitles = true;
		if (j == 48)
			config.use_lazy_dict = true;
		if (j == 49 || j == 53)
			config.store_query_quality = true;
	}
}

void print_staxids(TextBuffer &out, unsigned subject_global_id, const Metadata &metadata)
{
	out.print((*metadata.taxon_list)[subject_global_id], ';');
}

void Blast_tab_format::print_match(const Hsp_context& r, const Metadata &metadata, TextBuffer &out)
{
	for (vector<unsigned>::const_iterator i = fields.begin(); i != fields.end(); ++i) {
		switch (*i) {
		case 0:
			out.write_until(r.query_name, Const::id_delimiters);
			break;
		case 4:
			out << r.query.source().length();
			break;
		case 5:
			print_title(out, r.subject_name, false, false, "<>");
			break;
		case 6:
			print_title(out, r.subject_name, false, true, "<>");
			break;
		case 12:
			out << r.subject_len;
			break;
		case 13:
			out << r.oriented_query_range().begin_ + 1;
			break;
		case 14:
			out << r.oriented_query_range().end_ + 1;
			break;
		case 15:
			out << r.subject_range().begin_ + 1;
			break;
		case 16:
			out << r.subject_range().end_;
			break;
		case 17:
			r.query.source().print(out, r.query_source_range().begin_, r.query_source_range().end_, input_value_traits);
			break;
		case 18:
		{
			vector<Letter> seq;
			seq.reserve(r.subject_range().length());
			for (Hsp_context::Iterator j = r.begin(); j.good(); ++j)
				if (!(j.op() == op_insertion))
					seq.push_back(j.subject());
			out << sequence(seq);
			break;
		}
		case 19:
			out.print_e(r.evalue());
			break;
		case 20:
			out << r.bit_score();
			break;
		case 21:
			out << r.score();
			break;
		case 22:
			out << r.length();
			break;
		case 23:
			out << (double)r.identities() * 100 / r.length();
			break;
		case 24:
			out << r.identities();
			break;
		case 25:
			out << r.mismatches();
			break;
		case 26:
			out << r.positives();
			break;
		case 27:
			out << r.gap_openings();
			break;
		case 28:
			out << r.gaps();
			break;
		case 29:
			out << (double)r.positives() * 100.0 / r.length();
			break;
		case 31:
			out << r.blast_query_frame();
			break;
		case 33:
		{
			unsigned n_matches = 0;
			for (Hsp_context::Iterator i = r.begin(); i.good(); ++i) {
				switch (i.op()) {
				case op_match:
					++n_matches;
					break;
				case op_substitution:
				case op_frameshift_forward:
				case op_frameshift_reverse:
					if (n_matches > 0) {
						out << n_matches;
						n_matches = 0;
					}
					out << i.query_char() << i.subject_char();
					break;
				case op_insertion:
					if (n_matches > 0) {
						out << n_matches;
						n_matches = 0;
					}
					out << i.query_char() << '-';
					break;
				case op_deletion:
					if (n_matches > 0) {
						out << n_matches;
						n_matches = 0;
					}
					out << '-' << i.subject_char();
					break;
				}
			}
			if (n_matches > 0)
				out << n_matches;
		}
			break;
		case 34:
			print_staxids(out, r.orig_subject_id, metadata);
			break;
		case 35: {
			const vector<string> &names = *metadata.taxonomy_scientific_names;
			const vector<unsigned> tax_id = (*metadata.taxon_list)[r.orig_subject_id];
			for (size_t i = 0; i < tax_id.size(); ++i) {
				if (i > 0)
					out << ';';
				if (tax_id[i] < names.size() && !names[tax_id[i]].empty())
					out << names[tax_id[i]];
				else
					out << tax_id[i];
			}
			break;
		}
		case 39:
			print_title(out, r.subject_name, true, false, "<>");
			break;
		case 40:
			print_title(out, r.subject_name, true, true, "<>");
			break;
		case 43:
			out << (double)r.query_source_range().length()*100.0 / r.query.source().length();
			break;
		case 45:
			out << r.query_name;
			break;
		case 46:
			out << r.sw_score() - r.bit_score();
			break;
		case 47:
			out << r.time();
			break;
		case 48:
			out << r.subject_seq;
			break;
		case 49:
			out << (query_qual && (*query_qual)[r.query_id].present() ? (*query_qual)[r.query_id].substr(r.query_source_range().begin_, r.query_source_range().end_).c_str() : "*");
			break;
		case 50:
			out << query_block_to_database_id[r.query_id];
			break;
		case 51:
			out << r.orig_subject_id;
			break;
		case 52:
			out << (double)r.subject_range().length() * 100.0 / r.subject_len;
			break;
		case 53:
			out << (query_qual && (*query_qual)[r.query_id].present() ? (*query_qual)[r.query_id].c_str() : "*");
			break;
		case 54:
			r.query.source().print(out, input_value_traits);
			break;
		default:
			throw std::runtime_error(string("Invalid output field: ") + field_str[*i]);
		}
		if (i < fields.end() - 1)
			out << '\t';
	}
	out << '\n';
}

void Blast_tab_format::print_query_intro(size_t query_num, const char *query_name, unsigned query_len, TextBuffer &out, bool unaligned) const
{
	if (unaligned && config.report_unaligned == 1) {
		for (vector<unsigned>::const_iterator i = fields.begin(); i != fields.end(); ++i) {
			switch (*i) {
			case 0:
				out.write_until(query_name, Const::id_delimiters);
				break;
			case 4:
				out << query_len;
				break;
			case 5:
			case 6:
			case 17:
			case 18:
			case 33:
			case 35:
			case 39:
			case 40:
			case 48:
			case 49:
				out << '*';
				break;
			case 12:
			case 13:
			case 14:
			case 15:
			case 16:
			case 19:
			case 20:
			case 21:
			case 22:
			case 23:
			case 24:
			case 25:
			case 26:
			case 27:
			case 28:
			case 29:
			case 43:
				out << "-1";
				break;			
			case 31:
			case 34:
				out << '0';
				break;
			case 45:
				out << query_name;
				break;
			case 53:
				out << (query_qual && (*query_qual)[query_num].present() ? (*query_qual)[query_num].c_str() : "*");
				break;
			case 54:
				(align_mode.query_translated ? query_source_seqs::get()[query_num] : query_seqs::get()[query_num]).print(out, input_value_traits);
				break;
			default:
				throw std::runtime_error(string("Invalid output field: ") + field_str[*i]);
			}
			if (i < fields.end() - 1)
				out << '\t';			
		}
		out << '\n';
	}
}

void Blast_tab_format::print_header(Consumer &f, int mode, const char *matrix, int gap_open, int gap_extend, double evalue, const char *first_query_name, unsigned first_query_len) const {
	if (!config.output_header)
		return;
	stringstream ss;
	ss << "# DIAMOND v" << Const::version_string << ". http://github.com/bbuchfink/diamond" << endl;
	ss << "# Invocation: " << config.invocation << endl;
	ss << "# Fields: " << join(", ", apply(fields, [](unsigned i) -> string { return string(field_desc[i]); })) << endl;
	const string s(ss.str());
	f.consume(s.data(), s.length());
}
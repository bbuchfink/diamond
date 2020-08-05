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

#include <sstream>
#include <set>
#include <numeric>
#include <string.h>
#include <algorithm>
#include "../basic/match.h"
#include "output_format.h"
#include "../data/taxonomy.h"
#include "../data/queries.h"

using namespace std;

const vector<OutputField> Blast_tab_format::field_def = {
{ "qseqid", "Query ID", false, false },		// 0 means Query Seq - id
{ "qgi", "qgi", false, false },			// 1 means Query GI
{ "qacc", "qacc", false, false },			// 2 means Query accesion
{ "qaccver", "qaccver", false, false },		// 3 means Query accesion.version
{ "qlen", "Query length", false, false },			// 4 means Query sequence length
{ "sseqid",	"Subject ID", false, false },	// 5 means Subject Seq - id
{ "sallseqid", "Subject IDs", false, false },	// 6 means All subject Seq - id(s), separated by a ';'
{ "sgi", "sgi", false, false },			// 7 means Subject GI
{ "sallgi", "sallgi", false, false },		// 8 means All subject GIs
{ "sacc", "sacc", false, false },			// 9 means Subject accession
{ "saccver", "saccver", false, false },		// 10 means Subject accession.version
{ "sallacc", "sallacc", false, false },		// 11 means All subject accessions
{ "slen", "Subject length", false, false },			// 12 means Subject sequence length
{ "qstart", "Start of alignment in query", false, true },		// 13 means Start of alignment in query
{ "qend", "End of alignment in query", false, true },			// 14 means End of alignment in query
{ "sstart", "Start of alignment in subject", false, true },		// 15 means Start of alignment in subject
{ "send", "End of alignment in subject", false, true },			// 16 means End of alignment in subject
{ "qseq", "Aligned part of query sequence", false, true },			// 17 means Aligned part of query sequence
{ "sseq", "Aligned part of subject sequence", true, false },			// 18 means Aligned part of subject sequence
{ "evalue", "Expected value", false, false },		// 19 means Expect value
{ "bitscore", "Bit score", false, false },		// 20 means Bit score
{ "score", "Raw score", false, false },		// 21 means Raw score
{ "length", "Alignment length", false, true },		// 22 means Alignment length
{ "pident", "Percentage of identical matches", false, true },		// 23 means Percentage of identical matches
{ "nident", "Number of identical matches", false, true },		// 24 means Number of identical matches
{ "mismatch", "Number of mismatches", false, true },		// 25 means Number of mismatches
{ "positive", "Number of positive-scoring matches", false, true },		// 26 means Number of positive - scoring matches
{ "gapopen", "Number of gap openings", false, true },		// 27 means Number of gap openings
{ "gaps", "Total number of gaps", false, true },			// 28 means Total number of gaps
{ "ppos", "Percentage of positive-scoring matches", false, true },			// 29 means Percentage of positive - scoring matches
{ "frames", "frames", false, false },		// 30 means Query and subject frames separated by a '/'
{ "qframe", "Query frame", false, false },		// 31 means Query frame
{ "sframe", "sframe", false, false },		// 32 means Subject frame
{ "btop", "Blast traceback operations", true, false },			// 33 means Blast traceback operations(BTOP)
{ "staxids", "Subject Taxonomy IDs", false, false },		// 34 means unique Subject Taxonomy ID(s), separated by a ';'	(in numerical order)
{ "sscinames", "Subject scientific names", false, false },	// 35 means unique Subject Scientific Name(s), separated by a ';'
{ "scomnames", "scomnames", false, false },	// 36 means unique Subject Common Name(s), separated by a ';'
{ "sblastnames", "sblastnames", false, false },	// 37 means unique Subject Blast Name(s), separated by a ';'	(in alphabetical order)
{ "sskingdoms",	"Subject super kingdoms", false, false }, // 38 means unique Subject Super Kingdom(s), separated by a ';'	(in alphabetical order)
{ "stitle", "Subject title", false, false },		// 39 means Subject Title
{ "salltitles", "Subject titles", false, false },	// 40 means All Subject Title(s), separated by a '<>'
{ "sstrand", "sstrand", false, false },		// 41 means Subject Strand
{ "qcovs", "qcovs", false, true },		// 42 means Query Coverage Per Subject
{ "qcovhsp", "Query coverage per HSP", false, true },		// 43 means Query Coverage Per HSP
{ "qcovus", "qcovus", false, true },		// 44 means Query Coverage Per Unique Subject(blastn only)
{ "qtitle", "Query title", false, false },		// 45 means Query title
{ "swdiff", "swdiff", false, false },		// 46
{ "time", "time", false, false }, 		// 47
{ "full_sseq", "Subject sequence", false, false },	// 48
{ "qqual", "Aligned part of query quality values", false, true },		// 49
{ "qnum", "qnum", false, false },			// 50
{ "snum", "snum", false, false },			// 51
{ "scovhsp", "Subject coverage per HSP", false, true },		// 52
{ "full_qqual", "Query quality values", false, false},	// 53
{ "full_qseq", "Query sequence", false, false },	// 54
{ "qseq_gapped", "Query sequence with gaps", true, false },  // 55
{ "sseq_gapped", "Subject sequence with gaps", true, false },	// 56
{ "qstrand", "Query strand", false, false },		// 57
{ "cigar", "CIGAR", true, false },		// 58
{ "skingdoms", "Subject kingdoms", false, false },	// 59
{ "sphylums", "Subject phylums", false, false },		// 60
{ "ungapped_score", "Ungapped score", false, false }	// 61
};

Blast_tab_format::Blast_tab_format() :
	Output_format(blast_tab)
{
	static const unsigned stdf[] = { 0, 5, 23, 22, 25, 27, 13, 14, 15, 16, 19, 20 };
	const vector<string> &f = config.output_format;
	if (f.size() <= 1) {
		fields = vector<unsigned>(stdf, stdf + 12);
		if (config.frame_shift == 0) {
			needs_transcript = false;
			needs_stats = true;
		}
		return;
	}
	needs_transcript = false;
	needs_stats = false;
	for (vector<string>::const_iterator i = f.begin() + 1; i != f.end(); ++i) {
		auto it = std::find_if(field_def.begin(), field_def.end(), [i](const OutputField &f) {return f.key == *i; });
		if(it == field_def.end())
			throw std::runtime_error(string("Invalid output field: ") + *i);
		const auto j = it - field_def.begin();
		if (j == 34)
			needs_taxon_id_lists = true;
		if (j == 35 || j == 38 || j == 59 || j == 60) {
			needs_taxon_scientific_names = true;
			needs_taxon_id_lists = true;
		}
		if (j == 38 || j == 59 || j == 60) {
			needs_taxon_nodes = true;
			needs_taxon_ranks = true;
		}
		fields.push_back(j);
		if (j == 6 || j == 39 || j == 40 || j == 34)
			config.salltitles = true;
		if (j == 48)
			config.use_lazy_dict = true;
		if (j == 49 || j == 53)
			config.store_query_quality = true;
		if (field_def[j].need_transcript)
			needs_transcript = true;
		if (field_def[j].need_stats)
			needs_stats = true;
	}
	if (config.frame_shift != 0 && needs_stats)
		needs_transcript = true;
	if (config.traceback_mode == TracebackMode::NONE && config.max_hsps == 1 && !needs_transcript && !needs_stats && !config.query_range_culling && config.min_id == 0.0 && config.query_cover == 0.0 && config.subject_cover == 0.0)
		config.traceback_mode = TracebackMode::SCORE_ONLY;
}

void print_staxids(TextBuffer &out, unsigned subject_global_id, const Metadata &metadata)
{
	out.print((*metadata.taxon_list)[subject_global_id], ';');
}

template<typename _it>
void print_taxon_names(_it begin, _it end, const Metadata &metadata, TextBuffer &out) {
	if (begin == end) {
		out << "N/A";
		return;
	}
	const vector<string> &names = *metadata.taxonomy_scientific_names;
	for (_it i = begin; i != end; ++i) {
		if (i != begin)
			out << ';';
		if (*i < names.size() && !names[*i].empty())
			out << names[*i];
		else
			out << *i;
	}
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
			const vector<unsigned> tax_id = (*metadata.taxon_list)[r.orig_subject_id];
			print_taxon_names(tax_id.begin(), tax_id.end(), metadata, out);
			break;
		}
		case 38: {
			const set<unsigned> tax_id = metadata.taxon_nodes->rank_taxid((*metadata.taxon_list)[r.orig_subject_id], Rank::superkingdom);
			print_taxon_names(tax_id.begin(), tax_id.end(), metadata, out);
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
		case 49: {
			if (!query_qual) {
				out << '*';
				break;
			}
			const char *q = (*query_qual)[r.query_id];
			if (strlen(q) == 0) {
				out << '*';
				break;
			}
			out << string(q + r.query_source_range().begin_, q + r.query_source_range().end_).c_str();
			break;
		}			
		case 50:
			//out << (blocked_processing ? query_block_to_database_id[r.query_id] : r.query_id);
			out << r.query_id;
			break;
		case 51:
			out << r.orig_subject_id;
			break;
		case 52:
			out << (double)r.subject_range().length() * 100.0 / r.subject_len;
			break;
		case 53:
			out << (query_qual && strlen((*query_qual)[r.query_id]) ? (*query_qual)[r.query_id] : "*");
			break;
		case 54:
			r.query.source().print(out, input_value_traits);
			break;
		case 55:
			for (Hsp_context::Iterator i = r.begin(); i.good(); ++i)
				out << i.query_char();
			break;
		case 56:
			for (Hsp_context::Iterator i = r.begin(); i.good(); ++i)
				out << i.subject_char();
			break;
		case 57:
			if (align_mode.query_translated)
				out << ((r.blast_query_frame() > 0) ? '+' : '-');
			else
				out << '+';
			break;
		case 58:
			print_cigar(r, out);
			break;
		case 59: {
			const set<unsigned> tax_id = metadata.taxon_nodes->rank_taxid((*metadata.taxon_list)[r.orig_subject_id], Rank::kingdom);
			print_taxon_names(tax_id.begin(), tax_id.end(), metadata, out);
			break;
		}
		case 60: {
			const set<unsigned> tax_id = metadata.taxon_nodes->rank_taxid((*metadata.taxon_list)[r.orig_subject_id], Rank::phylum);
			print_taxon_names(tax_id.begin(), tax_id.end(), metadata, out);
			break;
		}
		case 61:
			out << score_matrix.bitscore(r.ungapped_score);
			break;
		default:
			throw std::runtime_error(string("Invalid output field: ") + field_def[*i].key);
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
			case 38:
			case 39:
			case 40:
			case 48:
			case 49:
			case 55:
			case 56:
			case 57:
			case 59:
			case 60:
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
			case 52:
			case 61:
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
				out << (query_qual && strlen((*query_qual)[query_num]) ? (*query_qual)[query_num] : "*");
				break;
			case 54:
				(align_mode.query_translated ? query_source_seqs::get()[query_num] : query_seqs::get()[query_num]).print(out, input_value_traits);
				break;
			default:
				throw std::runtime_error(string("Invalid output field: ") + field_def[*i].key);
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
	ss << "# Fields: " << join(", ", apply(fields, [](unsigned i) -> string { return string(field_def[i].description); })) << endl;
	const string s(ss.str());
	f.consume(s.data(), s.length());
}
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

#include <sstream>
#include <set>
#include <numeric>
#include <string.h>
#include <algorithm>
#include "../basic/match.h"
#include "output_format.h"
#include "../data/taxonomy.h"
#include "../data/queries.h"
#include "../util/util.h"
#include "../run/config.h"
#include "../util/sequence/sequence.h"
#include "../dp/ungapped.h"
#include "../basic/reduction.h"
#include "../dp/needleman_wunsch.h"

using namespace Output;
using std::endl;
using std::set;
using std::string;
using std::transform;
using std::back_inserter;
using std::vector;
using std::runtime_error;

const vector<OutputField> Blast_tab_format::field_def = {
{ "qseqid", "cseqid", "Query ID", HspValues::NONE, Flags::IS_STRING },		// 0 means Query Seq - id
{ "qgi", "", "qgi", HspValues::NONE, Flags::NONE },			// 1 means Query GI
{ "qacc", "", "qacc", HspValues::NONE, Flags::IS_STRING },			// 2 means Query accesion
{ "qaccver", "", "qaccver", HspValues::NONE, Flags::IS_STRING },		// 3 means Query accesion.version
{ "qlen", "clen", "Query length", HspValues::NONE, Flags::NONE },			// 4 means Query sequence length
{ "sseqid",	"mseqid", "Subject ID", HspValues::NONE, Flags::IS_STRING },	// 5 means Subject Seq - id
{ "sallseqid", "", "Subject IDs", HspValues::NONE, Flags::ALL_SEQIDS |Flags::IS_ARRAY },	// 6 means All subject Seq - id(s), separated by a ';'
{ "sgi", "", "sgi", HspValues::NONE, Flags::NONE },			// 7 means Subject GI
{ "sallgi", "", "sallgi", HspValues::NONE, Flags::ALL_SEQIDS | Flags::IS_ARRAY },		// 8 means All subject GIs
{ "sacc", "", "sacc", HspValues::NONE, Flags::IS_STRING },			// 9 means Subject accession
{ "saccver", "", "saccver", HspValues::NONE, Flags::IS_STRING },		// 10 means Subject accession.version
{ "sallacc", "", "sallacc", HspValues::NONE, Flags::ALL_SEQIDS | Flags::IS_STRING },		// 11 means All subject accessions
{ "slen", "mlen", "Subject length", HspValues::NONE, Flags::NONE },			// 12 means Subject sequence length
{ "qstart", "cstart", "Start of alignment in query", HspValues::QUERY_START, Flags::NONE },		// 13 means Start of alignment in query
{ "qend", "cend", "End of alignment in query", HspValues::QUERY_END, Flags::NONE },			// 14 means End of alignment in query
{ "sstart", "mstart", "Start of alignment in subject", HspValues::TARGET_START, Flags::NONE },		// 15 means Start of alignment in subject
{ "send", "mend", "End of alignment in subject", HspValues::TARGET_END, Flags::NONE },			// 16 means End of alignment in subject
{ "qseq", "", "Aligned part of query sequence", HspValues::QUERY_COORDS, Flags::IS_STRING },			// 17 means Aligned part of query sequence
{ "sseq", "", "Aligned part of subject sequence", HspValues::TRANSCRIPT, Flags::IS_STRING },			// 18 means Aligned part of subject sequence
{ "evalue", "evalue", "Expected value", HspValues::NONE, Flags::NONE },		// 19 means Expect value
{ "bitscore", "bitscore", "Bit score", HspValues::NONE, Flags::NONE },		// 20 means Bit score
{ "score", "score", "Raw score", HspValues::NONE, Flags::NONE },		// 21 means Raw score
{ "length", "length", "Alignment length", HspValues::LENGTH, Flags::NONE },		// 22 means Alignment length
{ "pident", "pident" ,"Percentage of identical matches", HspValues::IDENT | HspValues::LENGTH, Flags::NONE },		// 23 means Percentage of identical matches
{ "nident", "nident", "Number of identical matches", HspValues::IDENT, Flags::NONE },		// 24 means Number of identical matches
{ "mismatch", "mismatch", "Number of mismatches", HspValues::MISMATCHES, Flags::NONE },		// 25 means Number of mismatches
{ "positive", "positive", "Number of positive-scoring matches", HspValues::TRANSCRIPT, Flags::NONE },		// 26 means Number of positive - scoring matches
{ "gapopen", "gapopen", "Number of gap openings", HspValues::GAP_OPENINGS, Flags::NONE },		// 27 means Number of gap openings
{ "gaps", "gaps", "Total number of gaps", HspValues::GAPS, Flags::NONE },			// 28 means Total number of gaps
{ "ppos", "ppos", "Percentage of positive-scoring matches", HspValues::TRANSCRIPT, Flags::NONE },			// 29 means Percentage of positive - scoring matches
{ "frames", "", "frames", HspValues::NONE, Flags::IS_STRING },		// 30 means Query and subject frames separated by a '/'
{ "qframe", "", "Query frame", HspValues::NONE, Flags::NONE },		// 31 means Query frame
{ "sframe", "", "sframe", HspValues::NONE, Flags::NONE },		// 32 means Subject frame
{ "btop", "", "Blast traceback operations", HspValues::TRANSCRIPT, Flags::IS_STRING },			// 33 means Blast traceback operations(BTOP)
{ "staxids", "", "Subject Taxonomy IDs", HspValues::NONE, Flags::IS_ARRAY },		// 34 means unique Subject Taxonomy ID(s), separated by a ';'	(in numerical order)
{ "sscinames", "", "Subject scientific names", HspValues::NONE, Flags::IS_ARRAY },	// 35 means unique Subject Scientific Name(s), separated by a ';'
{ "scomnames", "", "scomnames", HspValues::NONE, Flags::IS_ARRAY },	// 36 means unique Subject Common Name(s), separated by a ';'
{ "sblastnames", "", "sblastnames", HspValues::NONE, Flags::NONE },	// 37 means unique Subject Blast Name(s), separated by a ';'	(in alphabetical order)
{ "sskingdoms",	"", "Subject super kingdoms", HspValues::NONE, Flags::IS_ARRAY }, // 38 means unique Subject Super Kingdom(s), separated by a ';'	(in alphabetical order)
{ "stitle", "", "Subject title", HspValues::NONE, Flags::FULL_TITLES | Flags::IS_STRING },		// 39 means Subject Title
{ "salltitles", "", "Subject titles", HspValues::NONE, Flags::FULL_TITLES | Flags::IS_ARRAY},	// 40 means All Subject Title(s), separated by a '<>'
{ "sstrand", "", "sstrand", HspValues::NONE, Flags::IS_STRING },		// 41 means Subject Strand
{ "qcovs", "", "qcovs", HspValues::NONE, Flags::NONE },		// 42 means Query Coverage Per Subject
{ "qcovhsp", "ccovhsp", "Query coverage per HSP", HspValues::QUERY_COORDS, Flags::NONE },		// 43 means Query Coverage Per HSP
{ "qcovus", "", "qcovus", HspValues::NONE, Flags::NONE },		// 44 means Query Coverage Per Unique Subject(blastn only)
{ "qtitle", "", "Query title", HspValues::NONE, Flags::IS_STRING },		// 45 means Query title
{ "swdiff", "", "swdiff", HspValues::NONE, Flags::NONE},		// 46
{ "time", "", "time", HspValues::NONE, Flags::NONE }, 		// 47
{ "full_sseq", "", "Subject sequence", HspValues::NONE, Flags::IS_STRING },	// 48
{ "qqual", "", "Aligned part of query quality values", HspValues::QUERY_COORDS, Flags::IS_STRING},		// 49
{ "qnum", "", "qnum", HspValues::NONE, Flags::NONE },			// 50
{ "snum", "", "snum", HspValues::NONE, Flags::NONE },			// 51
{ "scovhsp", "mcovhsp", "Subject coverage per HSP", HspValues::TARGET_COORDS, Flags::NONE },		// 52
{ "full_qqual", "", "Query quality values", HspValues::NONE, Flags::NONE },	// 53
{ "full_qseq", "", "Query sequence", HspValues::NONE, Flags::IS_STRING },	// 54
{ "qseq_gapped", "", "Query sequence with gaps", HspValues::TRANSCRIPT, Flags::IS_STRING },  // 55
{ "sseq_gapped", "", "Subject sequence with gaps", HspValues::TRANSCRIPT, Flags::IS_STRING },	// 56
{ "qstrand", "", "Query strand", HspValues::NONE, Flags::IS_STRING },		// 57
{ "cigar", "", "CIGAR", HspValues::TRANSCRIPT, Flags::IS_STRING },		// 58
{ "skingdoms", "", "Subject kingdoms", HspValues::NONE, Flags::IS_ARRAY },	// 59
{ "sphylums", "", "Subject phylums", HspValues::NONE, Flags::IS_ARRAY },		// 60
{ "ungapped_score", "", "Ungapped score", HspValues::NONE, Flags::NONE },	// 61
{ "full_qseq_mate", "", "Query sequence of the mate", HspValues::NONE, Flags::IS_STRING }, // 62
{ "qseq_translated", "", "Aligned part of query sequence (translated)", HspValues::TRANSCRIPT, Flags::IS_STRING }, // 63 needs transcript only in frameshift mode
{ "reduced_match_bitstring", "", "", HspValues::TRANSCRIPT, Flags::IS_STRING}, // 64
{ "normalized_bitscore_semiglobal", "", "", HspValues::NONE, Flags::NONE},  // 65
{ "hspnum", "", "", HspValues::NONE, Flags::NONE}, // 66
{ "normalized_bitscore_global", "", "", HspValues::NONE, Flags::SELF_ALN_SCORES}, // 67
{ "pident_global", "", "", HspValues::IDENT, Flags::NONE }, // 68
{ "ext_count", "", "", HspValues::NONE, Flags::NONE }, // 69
{ "nw_semiglobal", "", "", HspValues::NONE, Flags::NONE }, // 70
{ "approx_pident", "approx_pident", "", HspValues::COORDS, Flags::NONE },  // 71
{ "corrected_bitscore", "corrected_bitscore", "", HspValues::NONE, Flags::NONE }, // 72
{ "neg_evalue", "neg_evalue", "", HspValues::NONE, Flags::NONE }, // 73
{ "reserved1", "reserved1", "", HspValues::NONE, Flags::NONE}, // 74
{ "reserved2", "reserved2", "", HspValues::NONE, Flags::NONE} // 75
};

template<>
struct EnumTraits<Header> {
	static const SEMap<Header> from_string;
};

const SEMap<Header> EnumTraits<Header>::from_string = { {"0", Header::NONE}, {"simple", Header::SIMPLE}, {"verbose", Header::VERBOSE} };

Header Blast_tab_format::header_format(unsigned workflow) {
	const bool cluster = workflow == Config::cluster || workflow == Config::DEEPCLUST;
	if (workflow != Config::blastp && !cluster)
		throw runtime_error("header_format");
	if (!config.output_header.present())
		return Header::NONE;
	if (config.output_header.empty())
		return cluster ? Header::SIMPLE : Header::VERBOSE;
	if (config.output_header.size() > 1)
		throw runtime_error("Invalid header format: " + join(" ", config.output_header));
	const Header h = from_string<Header>(config.output_header.front());
	if (h == Header::VERBOSE && cluster)
		throw runtime_error("Verbose header format is not supported for cluster workflow.");
	return h;
}

Blast_tab_format::Blast_tab_format(bool json) :
	OutputFormat((json ? OutputFormat::json : OutputFormat::blast_tab), HspValues::NONE, Flags::NONE),
    is_json(json)
{
	static const unsigned stdf[] = { 0, 5, 23, 22, 25, 27, 13, 14, 15, 16, 19, 20 };
	const vector<string> &f = config.output_format;
	if (f.size() <= 1) {
		fields = vector<int64_t>(stdf, stdf + 12);
		if (config.frame_shift == 0)
			hsp_values = HspValues::QUERY_COORDS | HspValues::TARGET_COORDS | HspValues::LENGTH | HspValues::IDENT | HspValues::MISMATCHES | HspValues::GAP_OPENINGS;
		else
			hsp_values = HspValues::TRANSCRIPT;
		return;
	}
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
		if (j == 48 || j == 65 || j == 71)
			flags |= Flags::TARGET_SEQS;
		if (j == 49 || j == 53)
			config.store_query_quality = true;
		if (j == 62)
			needs_paired_end_info = true;
		if (j == 71 && score_matrix.name() != "BLOSUM62")
			throw std::runtime_error("Approximate identity is only supported for the BLOSUM62 scoring matrix.");
		if ((j == 62 || j == 63) && !align_mode.query_translated)
			throw std::runtime_error("Output field only supported for translated search.");
		hsp_values |= field_def[j].hsp_values;
		flags |= field_def[j].flags;
	}
	//if (config.traceback_mode == TracebackMode::NONE && config.max_hsps == 1 && !needs_transcript && !needs_stats && !config.query_range_culling && config.min_id == 0.0 && config.query_cover == 0.0 && config.subject_cover == 0.0)
		//config.traceback_mode = TracebackMode::SCORE_ONLY;
}

void print_staxids(TextBuffer &out, int64_t subject_global_id, const SequenceFile& db,bool json = false)
{
	out.print(db.taxids(subject_global_id), json ? ',':';');
}

template<typename _it>
void print_taxon_names(_it begin, _it end, const SequenceFile& db, TextBuffer &out, bool json = false) {
	if (begin == end) {
		out << "N/A";
		return;
	}
	for (_it i = begin; i != end; ++i) {
        if(json)
            out << "\"";
		if (i != begin)
			out << json ? ',': ';';
		out << db.taxon_scientific_name(*i);
        if(json)
            out << "\"";
	}
}

void Blast_tab_format::print_match(const HspContext& r, Output::Info& info)
{
	TextBuffer& out = info.out;
    const char* prepos = "\t";
    if(is_json) {
        if(r.hit_num !=0)
            out << ",";
        out << "\n\t{\n";
    }
	for (auto i = fields.cbegin(); i != fields.cend(); ++i) {
        if(is_json){
            out << prepos << "\"" << field_def[*i].key << "\":";
            if(flag_any(field_def[*i].flags , Flags::IS_STRING))
                out << "\"";
            if(flag_any(field_def[*i].flags , Flags::IS_ARRAY))
                out << "[";
        }
		switch (*i) {
		case 0:
			out.write_until(r.query_title.c_str(), Util::Seq::id_delimiters);
			break;
		case 4:
			out << r.query_len;
			break;
		case 5:
			print_title(out, r.target_title.c_str(), false, false, "<>");
			break;
		case 6:
            print_title(out, r.target_title.c_str(), false, true, is_json ? "," : "<>", 0,is_json);
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
			for (HspContext::Iterator j = r.begin(); j.good(); ++j)
				if (!(j.op() == op_insertion))
					seq.push_back(j.subject());
			out << Sequence(seq);
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
			out << r.id_percent();
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
			for (HspContext::Iterator i = r.begin(); i.good(); ++i) {
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
			print_staxids(out, r.subject_oid, *info.db, is_json);
			break;
		case 35: {
			const vector<TaxId> tax_id = info.db->taxids(r.subject_oid);
			print_taxon_names(tax_id.begin(), tax_id.end(), *info.db, out, is_json);
			break;
		}
		case 38: {
			const set<TaxId> tax_id = info.db->taxon_nodes().rank_taxid(info.db->taxids(r.subject_oid), Rank::superkingdom);
			print_taxon_names(tax_id.begin(), tax_id.end(), *info.db, out, is_json);
			break;
		}
		case 39:
            print_title(out, r.target_title.c_str(), true, false, is_json ? "," : "<>",0,is_json);
            break;
		case 40:
            print_title(out, r.target_title.c_str(), true, true, is_json ? ",":"<>", 0 ,is_json);
            break;
		case 43:
			out << r.qcovhsp();
			break;
		case 45:
                out << r.query_title;
			break;
		case 46:
			out << 0;
			break;
		case 47:
			out << 0;
			break;
		case 48:
			out << r.subject_seq;
			break;
		case 49: {
			if (strlen(info.query.qual) == 0) {
				out << '*';
				break;
			}
			out << string(info.query.qual + r.query_source_range().begin_, info.query.qual + r.query_source_range().end_).c_str();
			break;
		}
		case 50:
			out << r.query_oid;
			break;
		case 51:
			out << r.subject_oid;
			break;
		case 52:
			out << r.scovhsp();
			break;
		case 53:
			out << strlen(info.query.qual) ? info.query.qual : "*";
			break;
		case 54:
			r.query.source().print(out, input_value_traits);
			break;
		case 55:
			for (HspContext::Iterator i = r.begin(); i.good(); ++i)
				out << i.query_char();
			break;
		case 56:
			for (HspContext::Iterator i = r.begin(); i.good(); ++i)
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
			const set<TaxId> tax_id = info.db->taxon_nodes().rank_taxid(info.db->taxids(r.subject_oid), Rank::kingdom);
			print_taxon_names(tax_id.begin(), tax_id.end(), *info.db, out);
			break;
		}
		case 60: {
			const set<TaxId> tax_id = info.db->taxon_nodes().rank_taxid(info.db->taxids(r.subject_oid), Rank::phylum);
			print_taxon_names(tax_id.begin(), tax_id.end(), *info.db, out);
			break;
		}
		case 61:
			out << r.ungapped_score;
			break;
		case 62: {
			if (config.query_file.size() == 2) {
				info.query.mate_seq.print(out, input_value_traits);
				break;
			}
			else {
				out << '*';
				break;
			}
		}
		case 63: {
			if (config.frame_shift) {
				vector<Letter> seq;
				seq.reserve(r.query_range().length());
				for (HspContext::Iterator j = r.begin(); j.good(); ++j)
					if (j.op() != op_deletion && j.op() != op_frameshift_forward && j.op() != op_frameshift_reverse)
						seq.push_back(j.query());
				out << Sequence(seq);
			}
			else {
				r.query.index(r.frame()).print(out, r.query_range().begin_, r.query_range().end_, amino_acid_traits);
			}
			break;
		}
		case 64: {
			string s;
			size_t n = 0;
			for (HspContext::Iterator j = r.begin(); j.good(); ++j) {
				if (j.op() == op_deletion || j.op() == op_insertion) {
					if (!s.empty()) {
						if (n++ > 0) out << '\t';
						out << s;
						s.clear();
					}
				}
				else
					if (j.query() < 20 && j.subject() < 20) {
						if (Reduction::reduction(j.query()) == Reduction::reduction(j.subject()))
							s += '1';
						else
							s += '0';
					}
					else
						s += '0';
			}
			if (n > 0) out << '\t';
			out << s;
			break;
		}
		case 71:
			out << r.approx_id();
			break;
		case 72:
			out << r.corrected_bit_score();
			break;
#ifdef EXTRA
		case 65:
			out.print_d(r.bit_score() / score_matrix.bitscore(self_score(r.subject_seq)));
			break;
		case 66:
			out << r.hsp_num;
			break;
		case 67:
			out << r.bit_score() / std::max(r.query_self_aln_score, r.target_self_aln_score) * 100;
			break;
		case 68:
			out << (double)r.identities() / std::max(r.query.index(r.frame()).length(), r.subject_len) * 100;
			break;
		case 69:
			out << info.stats.extension_count;
			break;
		case 70:
			out.print_e(score_matrix.bitscore(r.query[Frame(0)].length() < r.subject_len ?
				nw_semiglobal(r.query[Frame(0)], r.subject_seq)
				: nw_semiglobal(r.subject_seq, r.query[Frame(0)])));
			break;
		case 73:
			out.print_e(r.evalue() == 0 ? r.bit_score() : -r.evalue());
			break;
		case 74:
			out << r.reserved1();
			break;
		case 75:
			out << r.reserved2();
			break;
#endif
		default:
			throw std::runtime_error(string("Invalid output field: ") + field_def[*i].key);
		}
        if(is_json) {
            if(flag_any(field_def[*i].flags , Flags::IS_STRING))
                out << "\"";
            if(flag_any(field_def[*i].flags , Flags::IS_ARRAY))
                out << "]";
            if ( i < fields.end() - 1)
                out << ",\n";
            else
                out << "\n";
        }
        else if (i < fields.end() - 1)
            out << '\t';
	}
    out << (is_json ? "\t}" : "\n");
}

void Blast_tab_format::print_query_intro(Output::Info& info) const
{
	TextBuffer& out = info.out;
	if (info.unaligned && config.report_unaligned == 1) {
		for (auto i = fields.cbegin(); i != fields.cend(); ++i) {
			switch (*i) {
			case 0:
				out.write_until(info.query.title, Util::Seq::id_delimiters);
				break;
			case 4:
				out << info.query.len;
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
			case 58:
			case 59:
			case 60:
			case 63:
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
			case 71:
				out << "-1";
				break;			
			case 31:
			case 34:
				out << '0';
				break;
			case 45:
				out << info.query.title;
				break;
			case 53:
				out << strlen(info.query.qual) ? info.query.qual : "*";
				break;
			case 54:
				info.query.source_seq.print(out, input_value_traits);
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

void Blast_tab_format::output_header(Consumer& f, bool cluster) const {
	vector<string> headers;
	transform(fields.begin(), fields.end(), back_inserter(headers), [cluster](int64_t i) {
		const string& key = cluster ? Blast_tab_format::field_def[i].clust_key : Blast_tab_format::field_def[i].key;
		if (cluster && key.empty())
			throw runtime_error("Output field not supported for clustering: " + Blast_tab_format::field_def[i].key);
		return key;
		});
	const string s = join("\t", headers) + '\n';
	f.consume(s.data(), s.length());
	return;
}

void Blast_tab_format::print_header(Consumer &f, int mode, const char *matrix, int gap_open, int gap_extend, double evalue, const char *first_query_name, unsigned first_query_len) const {
	const Header h = header_format(Config::blastp);
	if (h == Header::VERBOSE) {
		std::stringstream ss;
		ss << "# DIAMOND v" << Const::version_string << ". http://github.com/bbuchfink/diamond" << endl;
		ss << "# Invocation: " << config.invocation << endl;
		ss << "# Fields: " << join(", ", apply(fields, [](int64_t i) -> string { return string(field_def[i].description); })) << endl;
		const string s(ss.str());
		f.consume(s.data(), s.length());
	}
	else if (h == Header::SIMPLE) {
		output_header(f, false);
	}
    if(is_json) {
		const char* hdr = "[";
		f.consume(hdr, strlen(hdr));
    }
}
void Blast_tab_format::print_footer(Consumer &f) const {
    if(is_json) {
        const char* s = "\n]";
		f.consume(s, strlen(s));
    }
}
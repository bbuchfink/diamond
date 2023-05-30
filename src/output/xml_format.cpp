#include <sstream>
#include "output_format.h"
#include "../util/escape_sequences.h"
#include "../data/taxonomy.h"
#include "../util/util.h"
#include "../util/sequence/sequence.h"

using std::endl;
using std::string;

void XML_format::print_match(const HspContext& r, Output::Info& info)
{
	auto& out = info.out;
	if (r.hsp_num == 0) {
		if (r.hit_num > 0)
			out << "  </Hit_hsps>" << '\n' << "</Hit>" << '\n';

		out << "<Hit>" << '\n'
			<< "  <Hit_num>" << r.hit_num + 1 << "</Hit_num>" << '\n';
		string id, def;
		const string target_seqid = r.target_title;
		Util::Seq::get_title_def(target_seqid, id, def);
		if (config.xml_blord_format) {
			out << "  <Hit_id>gnl|BL_ORD_ID|" << r.subject_oid << "</Hit_id>" << '\n'
				<< "  <Hit_def>";
			OutputFormat::print_title(out, target_seqid.c_str(), true, true, " &gt;", &EscapeSequences::XML);
			out << "</Hit_def>" << '\n';
		}
		else {
			out << "  <Hit_id>";
			print_escaped(out, id, &EscapeSequences::XML);
			out << "</Hit_id>" << '\n'
				<< "  <Hit_def>";
			OutputFormat::print_title(out, def.c_str(), true, true, " &gt;", &EscapeSequences::XML);
			out << "</Hit_def>" << '\n';
		}
		out << "  <Hit_accession>";
		print_escaped(out, get_accession(id, info.acc_stats), &EscapeSequences::XML);
		out << "</Hit_accession>" << '\n'
			<< "  <Hit_len>" << r.subject_len << "</Hit_len>" << '\n'
			<< "  <Hit_hsps>" << '\n';
	}

	out << "    <Hsp>" << '\n'
		<< "      <Hsp_num>" << r.hsp_num + 1 << "</Hsp_num>" << '\n'
		<< "      <Hsp_bit-score>" << r.bit_score() << "</Hsp_bit-score>" << '\n'
		<< "      <Hsp_score>" << r.score() << "</Hsp_score>" << '\n'
		<< "      <Hsp_evalue>";
	out.print_e(r.evalue());
	out << "</Hsp_evalue>" << '\n'
		<< "      <Hsp_query-from>" << r.query_source_range().begin_ + 1 << "</Hsp_query-from>" << '\n'
		<< "      <Hsp_query-to>" << r.query_source_range().end_ << "</Hsp_query-to>" << '\n'
		<< "      <Hsp_hit-from>" << r.subject_range().begin_ + 1 << "</Hsp_hit-from>" << '\n'
		<< "      <Hsp_hit-to>" << r.subject_range().end_ << "</Hsp_hit-to>" << '\n'
		<< "      <Hsp_query-frame>" << r.blast_query_frame() << "</Hsp_query-frame>" << '\n'
		<< "      <Hsp_hit-frame>0</Hsp_hit-frame>" << '\n'
		<< "      <Hsp_identity>" << r.identities() << "</Hsp_identity>" << '\n'
		<< "      <Hsp_positive>" << r.positives() << "</Hsp_positive>" << '\n'
		<< "      <Hsp_gaps>" << r.gaps() << "</Hsp_gaps>" << '\n'
		<< "      <Hsp_align-len>" << r.length() << "</Hsp_align-len>" << '\n'
		<< "         <Hsp_qseq>";

	for (HspContext::Iterator i = r.begin(); i.good(); ++i)
		out << i.query_char();

	out << "</Hsp_qseq>" << '\n'
		<< "         <Hsp_hseq>";

	for (HspContext::Iterator i = r.begin(); i.good(); ++i)
		out << i.subject_char();

	out << "</Hsp_hseq>" << '\n'
		<< "      <Hsp_midline>";

	for (HspContext::Iterator i = r.begin(); i.good(); ++i)
		out << i.midline_char(score_matrix(i.query(), i.subject()));

	out << "</Hsp_midline>" << '\n'
		<< "    </Hsp>" << '\n';
}

void XML_format::print_header(Consumer& f, int mode, const char* matrix, int gap_open, int gap_extend, double evalue, const char* first_query_name, unsigned first_query_len) const
{
	std::stringstream ss;
	ss << "<?xml version=\"1.0\"?>" << endl
		<< "<!DOCTYPE BlastOutput PUBLIC \"-//NCBI//NCBI BlastOutput/EN\" \"http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd\">" << endl
		<< "<BlastOutput>" << endl
		<< "  <BlastOutput_program>" << mode_str(mode) << "</BlastOutput_program>" << endl
		<< "  <BlastOutput_version>" << Const::program_name << ' ' << Const::version_string << "</BlastOutput_version>" << endl
		<< "  <BlastOutput_reference>Benjamin Buchfink, Xie Chao, and Daniel Huson (2015), &quot;Fast and sensitive protein alignment using DIAMOND&quot;, Nature Methods 12:59-60.</BlastOutput_reference>" << endl
		<< "  <BlastOutput_db>" << config.database << "</BlastOutput_db>" << endl
		<< "  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>" << endl
		<< "  <BlastOutput_query-def>";
	const string fqn = string(first_query_name);
	string escaped;
	EscapeSequences::XML.escape(fqn, escaped);
	ss << escaped.substr(0, fqn.find('\1'));
	ss << "</BlastOutput_query-def>" << endl << "  <BlastOutput_query-len>" << first_query_len << "</BlastOutput_query-len>" << endl
		<< "  <BlastOutput_param>" << endl
		<< "    <Parameters>" << endl
		<< "      <Parameters_matrix>" << matrix << "</Parameters_matrix>" << endl
		<< "      <Parameters_expect>" << evalue << "</Parameters_expect>" << endl
		<< "      <Parameters_gap-open>" << gap_open << "</Parameters_gap-open>" << endl
		<< "      <Parameters_gap-extend>" << gap_extend << "</Parameters_gap-extend>" << endl
		<< "      <Parameters_filter>F</Parameters_filter>" << endl
		<< "    </Parameters>" << endl
		<< "  </BlastOutput_param>" << endl
		<< "<BlastOutput_iterations>" << endl;
	f.consume(ss.str().c_str(), ss.str().length());
}

void XML_format::print_query_intro(Output::Info& info) const
{
	info.out << "<Iteration>" << '\n'
		<< "  <Iteration_iter-num>" << info.query.oid + 1 << "</Iteration_iter-num>" << '\n'
		<< "  <Iteration_query-ID>Query_" << info.query.oid + 1 << "</Iteration_query-ID>" << '\n'
		<< "  <Iteration_query-def>";
	print_title(info.out, info.query.title, true, false, "", &EscapeSequences::XML);
	info.out << "</Iteration_query-def>" << '\n'
		<< "  <Iteration_query-len>" << info.query.len << "</Iteration_query-len>" << '\n'
		<< "<Iteration_hits>" << '\n';
}

void XML_format::print_query_epilog(Output::Info& info) const
{
	if (!info.unaligned) {
		info.out << "  </Hit_hsps>" << '\n'
			<< "</Hit>" << '\n';
	}
	((info.out << "</Iteration_hits>" << '\n'
		<< "  <Iteration_stat>" << '\n'
		<< "    <Statistics>" << '\n'
		<< "      <Statistics_db-num>" << info.db->sequence_count() << "</Statistics_db-num>" << '\n'
		<< "      <Statistics_db-len>" << info.db->letters() << "</Statistics_db-len>" << '\n'
		<< "      <Statistics_hsp-len>0</Statistics_hsp-len>" << '\n'
		<< "      <Statistics_eff-space>0</Statistics_eff-space>" << '\n'
		<< "      <Statistics_kappa>").print_d(score_matrix.k()) << "</Statistics_kappa>" << '\n'
		<< "      <Statistics_lambda>").print_d(score_matrix.lambda()) << "</Statistics_lambda>" << '\n'
		<< "      <Statistics_entropy>0</Statistics_entropy>" << '\n'
		<< "    </Statistics>" << '\n'
		<< "  </Iteration_stat>" << '\n'
		<< "</Iteration>" << '\n';
}

void XML_format::print_footer(Consumer& f) const
{
	std::stringstream ss;
	ss << "</BlastOutput_iterations>" << endl
		<< "</BlastOutput>";
	f.consume(ss.str().c_str(), ss.str().length());
}
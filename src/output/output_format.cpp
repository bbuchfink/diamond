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

#include <iostream>
#include "../data/taxonomy.h"
#include "output_format.h"
#include "../data/reference.h"

using std::endl;

auto_ptr<Output_format> output_format;

void Output_format::print_title(Text_buffer &buf, const char *id, bool full_titles, bool all_titles, const char *separator)
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
	const vector<string> t(tokenize(id, "\1"));
	vector<string>::const_iterator i = t.begin();
	for (; i<t.end() - 1; ++i) {
		if (full_titles)
			buf << *i << separator;
		else {
			buf.write_until(i->c_str(), Const::id_delimiters);
			buf << ";";
		}
		//n += i->length() + 2;
	}
	if (full_titles)
		buf << *i;
	else
		buf.write_until(i->c_str(), Const::id_delimiters);
	//n += i->length();
	//return n;
}

void print_hsp(Hsp_data &hsp, sequence query)
{
	Text_buffer buf;
	Pairwise_format().print_match(Hsp_context(hsp, 0, query, query, "", 0, 0, "", 0, 0, 0), buf);
	buf << '\0';
	cout << buf.get_begin();
}

Output_format* get_output_format()
{
	const vector<string> &f = config.output_format;
	if (f.size() == 0) {
		if (config.daa_file == "" || config.command == Config::view)
			return new Blast_tab_format;
		else if ((config.command == Config::blastp || config.command == Config::blastx) && config.daa_file.length() > 0)
			return new DAA_format();
	}
	if (f[0] == "tab" || f[0] == "6")
		return new Blast_tab_format;
	else if (f[0] == "sam" || f[0] == "101")
		return new Sam_format;
	else if (f[0] == "xml" || f[0] == "5")
		return new XML_format;
	else if ((config.command == Config::blastp || config.command == Config::blastx) && (f[0] == "daa" || f[0] == "100"))
		return new DAA_format;
	else if (f[0] == "0")
		return new Pairwise_format;
	else if (f[0] == "null")
		return new Null_format;
	else if (f[0] == "102")
		return new Taxon_format;
	else
		throw std::runtime_error("Invalid output format. Allowed values: 0,5,6,100,101,102");
}

void init_output()
{
	output_format = auto_ptr<Output_format>(get_output_format());
	if (*output_format == Output_format::taxon && config.toppercent == 100.0)
		config.toppercent = 10.0;
	if (config.toppercent == 100.0) {
		message_stream << "#Target sequences to report alignments for: ";
		if (config.max_alignments == 0) {
			config.max_alignments = std::numeric_limits<uint64_t>::max();
			message_stream << "unlimited" << endl;
		}
		else
			message_stream << config.max_alignments << endl;
	}
	else
		message_stream << "Percentage range of top alignment score to report hits: " << config.toppercent << endl;
}

void XML_format::print_match(const Hsp_context &r, Text_buffer &out)
{
	if(r.hsp_num == 0) {
		if (r.hit_num > 0)
			out << "  </Hit_hsps>" << '\n' << "</Hit>" << '\n';
		string id, def;
		get_title_def(r.subject_name, id, def);
		out << "<Hit>" << '\n'
			<< "  <Hit_num>" << r.hit_num + 1 << "</Hit_num>" << '\n'
			<< "  <Hit_id>" << id << "</Hit_id>" << '\n'
			<< "  <Hit_def>" << def << "</Hit_def>" << '\n';
		out << "  <Hit_accession>" << get_accession(id) << "</Hit_accession>" << '\n'
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
		<< "      <Hsp_query-to>" << r.query_source_range().end_ + 1 << "</Hsp_query-to>" << '\n'
		<< "      <Hsp_hit-from>" << r.subject_range().begin_ + 1 << "</Hsp_hit-from>" << '\n'
		<< "      <Hsp_hit-to>" << r.subject_range().end_ << "</Hsp_hit-to>" << '\n'
		<< "      <Hsp_query-frame>" << r.blast_query_frame() << "</Hsp_query-frame>" << '\n'
		<< "      <Hsp_hit-frame>0</Hsp_hit-frame>" << '\n'
		<< "      <Hsp_identity>" << r.identities() << "</Hsp_identity>" << '\n'
		<< "      <Hsp_positive>" << r.positives() << "</Hsp_positive>" << '\n'
		<< "      <Hsp_gaps>" << r.gaps() << "</Hsp_gaps>" << '\n'
		<< "      <Hsp_align-len>" << r.length() << "</Hsp_align-len>" << '\n'
		<< "         <Hsp_qseq>";

	for (Hsp_context::Iterator i = r.begin(); i.good(); ++i)
		out << i.query_char();
		
	out << "</Hsp_qseq>" << '\n'
		<< "         <Hsp_hseq>";

	for (Hsp_context::Iterator i = r.begin(); i.good(); ++i)
		out << i.subject_char();

	out << "</Hsp_hseq>" << '\n'
		<< "      <Hsp_midline>";

	for (Hsp_context::Iterator i = r.begin(); i.good(); ++i)
		out << i.midline_char();

	out << "</Hsp_midline>" << '\n'
		<< "    </Hsp>" << '\n';
}

void XML_format::print_header(Output_stream &f, int mode, const char *matrix, int gap_open, int gap_extend, double evalue, const char *first_query_name, unsigned first_query_len) const
{
	std::stringstream ss;
	ss << "<?xml version=\"1.0\"?>" << endl
		<< "<!DOCTYPE BlastOutput PUBLIC \"-//NCBI//NCBI BlastOutput/EN\" \"http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd\">" << endl
		<< "<BlastOutput>" << endl
		<< "  <BlastOutput_program>" << mode_str(mode) << "</BlastOutput_program>" << endl
		<< "  <BlastOutput_version>" << Const::program_name << ' ' << Const::version_string << "</BlastOutput_version>" << endl
		<< "  <BlastOutput_reference>Benjamin Buchfink, Xie Chao, and Daniel Huson (2015), &quot;Fast and sensitive protein alignment using DIAMOND&quot;, Nature Methods 12:59-60.</BlastOutput_reference>" << endl
		<< "  <BlastOutput_db></BlastOutput_db>" << endl
		<< "  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>" << endl
		<< "  <BlastOutput_query-def>";
	const string fqn = string(first_query_name);
	ss << fqn.substr(0, fqn.find('\1'));
		ss << "</BlastOutput_query-def>" << endl << "  <BlastOutput_query-len>" << first_query_len << "</BlastOutput_query-len>" << endl
		<< "  <BlastOutput_param>" << endl
		<< "    <Parameters>" << endl
		<< "      <Parameters_matrix>" << matrix << "</Parameters_matrix>" << endl
		<< "      <Parameters_expect>" << evalue << "</Parameters_expect>" << endl
		<< "      <Parameters_gap-open>" << gap_open << "</Parameters_gap-open>" << endl
		<< "      <Parameters_gap-extend>" << gap_extend << "</Parameters_gap-extend>" << endl
		<< "      <Parameters_filter></Parameters_filter>" << endl
		<< "    </Parameters>" << endl
		<< "  </BlastOutput_param>" << endl
		<< "<BlastOutput_iterations>" << endl;
	f.write(ss.str().c_str(), ss.str().length());
}

void XML_format::print_query_intro(size_t query_num, const char *query_name, unsigned query_len, Text_buffer &out, bool unaligned) const
{
	out << "<Iteration>" << '\n'
		<< "  <Iteration_iter-num>" << query_num + 1 << "</Iteration_iter-num>" << '\n'
		<< "  <Iteration_query-ID>Query_" << query_num + 1 << "</Iteration_query-ID>" << '\n'
		<< "  <Iteration_query-def>";
	print_title(out, query_name, true, false, "");
	out << "</Iteration_query-def>" << '\n'
		<< "  <Iteration_query-len>" << query_len << "</Iteration_query-len>" << '\n'
		<< "<Iteration_hits>" << '\n';
}

void XML_format::print_query_epilog(Text_buffer &out, const char *query_title, bool unaligned) const
{
	if (!unaligned) {
		out << "  </Hit_hsps>" << '\n'
			<< "</Hit>" << '\n';
	}
	((out << "</Iteration_hits>" << '\n'
		<< "  <Iteration_stat>" << '\n'
		<< "    <Statistics>" << '\n'
		<< "      <Statistics_db-num>" << (size_t)ref_header.sequences << "</Statistics_db-num>" << '\n'
		<< "      <Statistics_db-len>" << (size_t)config.db_size << "</Statistics_db-len>" << '\n'
		<< "      <Statistics_hsp-len>0</Statistics_hsp-len>" << '\n'
		<< "      <Statistics_eff-space>0</Statistics_eff-space>" << '\n'
		<< "      <Statistics_kappa>").print_d(score_matrix.k()) << "</Statistics_kappa>" << '\n'
		<< "      <Statistics_lambda>").print_d(score_matrix.lambda()) << "</Statistics_lambda>" << '\n'
		<< "      <Statistics_entropy>0</Statistics_entropy>" << '\n'
		<< "    </Statistics>" << '\n'
		<< "  </Iteration_stat>" << '\n'
		<< "</Iteration>" << '\n';
}

void XML_format::print_footer(Output_stream &f) const
{
	std::stringstream ss;
	ss << "</BlastOutput_iterations>" << endl
		<< "</BlastOutput>";
	f.write(ss.str().c_str(), ss.str().length());
}
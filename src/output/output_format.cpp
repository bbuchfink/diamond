/****
Copyright (c) 2016, Benjamin Buchfink
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

#include <iostream>
#include "output_format.h"

using std::endl;

void XML_format::print_match(const DAA_query_record::Match &r, Text_buffer &out) const
{
	if(r.hsp_num == 0) {
		if (r.hit_num > 0)
			out << "  </Hit_hsps>" << '\n' << "</Hit>" << '\n';
		out << "<Hit>" << '\n'
			<< "  <Hit_num>" << r.hit_num+1 << "</Hit_num>" << '\n'
			<< "  <Hit_id></Hit_id>" << '\n'
			<< "  <Hit_def>";
		if (r.subject_name.find_first_of('\1') == string::npos)
			out << r.subject_name;
		else
			this->print_salltitles(out, r.subject_name.c_str());
		out <<	"</Hit_def> " << '\n'
			<< "  <Hit_accession></Hit_accession>" << '\n'
			<< "  <Hit_len>" << r.total_subject_len << "</Hit_len>" << '\n'
			<< "  <Hit_hsps>" << '\n';
	}

	out << "    <Hsp>" << '\n'
		<< "      <Hsp_num>" << r.hsp_num + 1 << "</Hsp_num>" << '\n'
		<< "      <Hsp_bit-score>" << r.bit_score() << "</Hsp_bit-score>" << '\n'
		<< "      <Hsp_score>" << r.score << "</Hsp_score>" << '\n'
		<< "      <Hsp_evalue>";
	out.print_e(r.evalue());
	out << "</Hsp_evalue>" << '\n'
		<< "      <Hsp_query-from>" << r.unoriented_query_range().first + 1 << "</Hsp_query-from>" << '\n'
		<< "      <Hsp_query-to>" << r.unoriented_query_range().second + 1 << "</Hsp_query-to>" << '\n'
		<< "      <Hsp_hit-from>" << r.subject_begin + 1 << "</Hsp_hit-from>" << '\n'
		<< "      <Hsp_hit-to>" << r.subject_begin + r.subject_len << "</Hsp_hit-to>" << '\n'
		<< "      <Hsp_query-frame>" << r.blast_query_frame() << "</Hsp_query-frame>" << '\n'
		<< "      <Hsp_hit-frame>0</Hsp_hit-frame>" << '\n'
		<< "      <Hsp_identity>" << r.identities << "</Hsp_identity>" << '\n'
		<< "      <Hsp_positive>" << r.positives << "</Hsp_positive>" << '\n'
		<< "      <Hsp_gaps>" << r.gaps << "</Hsp_gaps>" << '\n'
		<< "      <Hsp_align-len>" << r.len << "</Hsp_align-len>" << '\n'
		<< "         <Hsp_qseq>";

	for (DAA_query_record::Match::Position_iterator i = r.begin(); i.good(); ++i)
		out << i.query_char();
		
	out << "</Hsp_qseq>" << '\n'
		<< "         <Hsp_hseq>";

	for (DAA_query_record::Match::Position_iterator i = r.begin(); i.good(); ++i)
		out << i.subject_char();

	out << "</Hsp_hseq>" << '\n'
		<< "      <Hsp_midline>";

	for (DAA_query_record::Match::Position_iterator i = r.begin(); i.good(); ++i)
		out << i.midline_char();

	out << "</Hsp_midline>" << '\n'
		<< "    </Hsp>" << '\n';
}

void XML_format::print_header(Output_stream &f, int mode, const char *matrix, int gap_open, int gap_extend, double evalue) const
{
	std::stringstream ss;
	ss << "<?xml version=\"1.0\"?>" << endl
		<< "<!DOCTYPE BlastOutput PUBLIC \"-//NCBI//NCBI BlastOutput/EN\" \"http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd\">" << endl
		<< "<BlastOutput>" << endl
		<< "  <BlastOutput_program>" << mode_str(mode) << "</BlastOutput_program>" << endl
		<< "  <BlastOutput_version>" << Const::program_name << ' ' << Const::version_string << "</BlastOutput_version>" << endl
		<< "  <BlastOutput_reference>Benjamin Buchfink, Xie Chao, and Daniel Huson (2015), &quot;Fast and sensitive protein alignment using DIAMOND&quot;, Nature Methods 12:59-60.</BlastOutput_reference>" << endl
		<< "  <BlastOutput_db></BlastOutput_db>" << endl
		<< "  <BlastOutput_query-ID></BlastOutput_query-ID>" << endl
		<< "  <BlastOutput_query-def></BlastOutput_query-def>" << endl
		<< "  <BlastOutput_query-len></BlastOutput_query-len>" << endl
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

void XML_format::print_query_intro(const DAA_query_record &r, Text_buffer &out) const
{
	out << "<Iteration>" << '\n'
		<< "  <Iteration_iter-num>" << r.query_num+1 << "</Iteration_iter-num>" << '\n'
		<< "  <Iteration_query-ID>Query_" << r.query_num+1 << "</Iteration_query-ID>" << '\n'
		<< "  <Iteration_query-def>" << r.query_name << "</Iteration_query-def>" << '\n'
		<< "  <Iteration_query-len>" << r.query_len() << "</Iteration_query-len>" << '\n'
		<< "<Iteration_hits>" << '\n';
}

void XML_format::print_query_epilog(const DAA_query_record &r, Text_buffer &out) const
{
	((out << "  </Hit_hsps>" << '\n'
		<< "</Hit>" << '\n'
		<< "</Iteration_hits>" << '\n'
		<< "  <Iteration_stat>" << '\n'
		<< "    <Statistics>" << '\n'
		<< "      <Statistics_db-num></Statistics_db-num>" << '\n'
		<< "      <Statistics_db-len>" << r.db_letters() << "</Statistics_db-len>" << '\n'
		<< "      <Statistics_hsp-len></Statistics_hsp-len>" << '\n'
		<< "      <Statistics_eff-space></Statistics_eff-space>" << '\n'
		<< "      <Statistics_kappa>").print_d(r.kappa()) << "</Statistics_kappa>" << '\n'
		<< "      <Statistics_lambda>").print_d(r.lambda()) << "</Statistics_lambda>" << '\n'
		<< "      <Statistics_entropy></Statistics_entropy>" << '\n'
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
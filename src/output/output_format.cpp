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

#include <iostream>
#include <limits.h>
#include <algorithm>
#include "../data/taxonomy.h"
#include "output_format.h"
#include "../data/reference.h"
#include "../util/escape_sequences.h"

using namespace std;

unique_ptr<Output_format> output_format;

void Output_format::print_title(TextBuffer &buf, const char *id, bool full_titles, bool all_titles, const char *separator, const EscapeSequences *esc)
{
	if (!all_titles) {
		print_escaped_until(buf, id, full_titles ? "\1" : Const::id_delimiters, esc);
		return;
	}
	if (strchr(id, '\1') == 0) {
		print_escaped_until(buf, id, full_titles ? "\1" : Const::id_delimiters, esc);
		return;
	}
	const vector<string> t(tokenize(id, "\1"));
	vector<string>::const_iterator i = t.begin();
	for (; i < t.end() - 1; ++i) {
		if (full_titles) {
			print_escaped(buf, *i, esc);
			buf << separator;
		}
		else {
			print_escaped_until(buf, i->c_str(), Const::id_delimiters, esc);
			buf << ";";
		}
	}
	if (full_titles)
		print_escaped(buf, *i, esc);
	else
		print_escaped_until(buf, i->c_str(), Const::id_delimiters, esc);
}

void print_hsp(Hsp &hsp, const TranslatedSequence &query)
{
	TextBuffer buf;
	Pairwise_format().print_match(Hsp_context(hsp, 0, query, "", 0, 0, "", 0, 0, 0, sequence()), Metadata(), buf);
	buf << '\0';
	cout << buf.get_begin() << endl;
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
	else if (f[0] == "daa" || f[0] == "100")
		return new DAA_format;
	else if (f[0] == "0")
		return new Pairwise_format;
	else if (f[0] == "null")
		return new Null_format;
	else if (f[0] == "102")
		return new Taxon_format;
	else if (f[0] == "paf" || f[0] == "103")
		return new PAF_format;
	else if (f[0] == "bin1")
		return new Bin1_format;
	else if (f[0] == "clus")
		return new Clustering_format(&f[1]);
	else
		throw std::runtime_error("Invalid output format. Allowed values: 0,5,6,100,101,102");
}

void init_output(bool have_taxon_id_lists, bool have_taxon_nodes, bool have_taxon_scientific_names)
{
	output_format = unique_ptr<Output_format>(get_output_format());
	if(config.command == Config::view && (output_format->needs_taxon_id_lists || output_format->needs_taxon_nodes || output_format->needs_taxon_scientific_names))
		throw runtime_error("Taxonomy features are not supported for the DAA format.");
	if (output_format->needs_taxon_id_lists && !have_taxon_id_lists)
		throw runtime_error("Output format requires taxonomy mapping information built into the database (use --taxonmap parameter for the makedb command).");
	if (output_format->needs_taxon_nodes && !have_taxon_nodes)
		throw runtime_error("Output format requires taxonomy nodes information built into the database (use --taxonnodes parameter for the makedb command).");
	if (output_format->needs_taxon_scientific_names && !have_taxon_scientific_names)
		throw runtime_error("Output format requires taxonomy names information built into the database (use --taxonnames parameter for the makedb command).");
	if (*output_format == Output_format::taxon && config.toppercent == 100.0 && config.min_bit_score == 0.0)
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
	if (config.traceback_mode == TracebackMode::NONE)
		config.traceback_mode = TracebackMode::VECTOR;
	log_stream << "Traceback mode: " << (int)config.traceback_mode << endl;
	log_stream << "Format options: transcript=" << output_format->needs_transcript << " stats=" << output_format->needs_stats << endl;
}

void Bin1_format::print_query_intro(size_t query_num, const char *query_name, unsigned query_len, TextBuffer &out, bool unaligned) const {
	out.write(std::numeric_limits<uint32_t>::max());
	out.write((uint32_t)query_num);
}

void Bin1_format::print_match(const Hsp_context& r, const Metadata &metadata, TextBuffer &out) {
	if (r.query_id < r.subject_id) {
		out.write((uint32_t)r.subject_id);
		out.write(r.bit_score() / std::max((unsigned)r.query.source().length(), r.subject_len));
	}
}


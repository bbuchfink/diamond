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
#include "../data/queries.h"
#include "output.h"
#include "../util/util.h"
#include "../run/config.h"
#include "../util/sequence/sequence.h"

using namespace std;

unique_ptr<Output_format> output_format;

void IntermediateRecord::read(BinaryBuffer::Iterator& f)
{
	f.read(target_dict_id);
	/*if (config.global_ranking_targets > 0) {
		uint16_t s;
		f.read(s);
		score = s;
		return;
	}*/
	if (*output_format == Output_format::daa)
		f.read(target_oid);
	f.read(flag);
	f.read_packed(flag & 3, score);
	f.read(evalue);

	if (output_format->hsp_values == Output::NONE)
		return;

	f.read_packed((flag >> 2) & 3, query_begin);
	f.read_varint(query_end);
	f.read_packed((flag >> 4) & 3, subject_begin);

	if (output_format->hsp_values & Output::TRANSCRIPT)
		transcript.read(f);
	else {
		f.read_varint(subject_end);
		f.read_varint(identities);
		f.read_varint(mismatches);
		f.read_varint(positives);
		f.read_varint(gap_openings);
		f.read_varint(gaps);
	}
}

interval IntermediateRecord::absolute_query_range() const
{
	if (query_begin < query_end)
		return interval(query_begin, query_end + 1);
	else
		return interval(query_end, query_begin + 1);
}

size_t IntermediateRecord::write_query_intro(TextBuffer& buf, unsigned query_id)
{
	size_t seek_pos = buf.size();
	buf.write((uint32_t)query_id).write((uint32_t)0);
	return seek_pos;
}

void IntermediateRecord::finish_query(TextBuffer& buf, size_t seek_pos)
{
	*(uint32_t*)(&buf[seek_pos + sizeof(uint32_t)]) = safe_cast<uint32_t>(buf.size() - seek_pos - sizeof(uint32_t) * 2);
}

void IntermediateRecord::write(TextBuffer& buf, const Hsp& match, unsigned query_id, size_t target_dict_id, size_t target_oid)
{
	const interval oriented_range(match.oriented_range());
	buf.write((uint32_t)target_dict_id);
	if (*output_format == Output_format::daa)
		buf.write((uint32_t)target_oid);
	buf.write(get_segment_flag(match));
	buf.write_packed(match.score);
	buf.write(match.evalue);
	if (output_format->hsp_values == Output::NONE)
		return;

	buf.write_packed(oriented_range.begin_);
	buf.write_varint(oriented_range.end_);
	buf.write_packed(match.subject_range.begin_);

	if (output_format->hsp_values & Output::TRANSCRIPT)
		buf << match.transcript.data();
	else {
		buf.write_varint(match.subject_range.end_);
		buf.write_varint(match.identities);
		buf.write_varint(match.mismatches);
		buf.write_varint(match.positives);
		buf.write_varint(match.gap_openings);
		buf.write_varint(match.gaps);
	}
}

void IntermediateRecord::write(TextBuffer& buf, uint32_t target_block_id, int score, const Search::Config& cfg) {
	const uint32_t target_oid = (uint32_t)cfg.target->block_id2oid(target_block_id);
	assert(target_oid < cfg.db_seqs);
	buf.write(target_oid);
	const uint16_t s = (uint16_t)std::min(score, USHRT_MAX);
	buf.write(s);
}

void IntermediateRecord::finish_file(Consumer& f)
{
	const uint32_t i = FINISHED;
	f.consume(reinterpret_cast<const char*>(&i), 4);
}

void Output_format::print_title(TextBuffer &buf, const char *id, bool full_titles, bool all_titles, const char *separator, const EscapeSequences *esc)
{
	if (!all_titles) {
		if (config.short_seqids)
			buf << Util::Seq::seqid(id, true);
		else
			print_escaped_until(buf, id, full_titles ? "\1" : Util::Seq::id_delimiters, esc);
		return;
	}
	if (strchr(id, '\1') == 0) {
		print_escaped_until(buf, id, full_titles ? "\1" : Util::Seq::id_delimiters, esc);
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
			print_escaped_until(buf, i->c_str(), Util::Seq::id_delimiters, esc);
			buf << ";";
		}
	}
	if (full_titles)
		print_escaped(buf, *i, esc);
	else
		print_escaped_until(buf, i->c_str(), Util::Seq::id_delimiters, esc);
}

void print_hsp(Hsp &hsp, const TranslatedSequence &query)
{
	TextBuffer buf;
	//Pairwise_format().print_match(Hsp_context(hsp, 0, query, "", 0, 0, 0, 0, Sequence()), Search::Config(true), buf);
	buf << '\0';
	cout << buf.data() << endl;
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
	else if (f[0] == "bin")
		return new Binary_format;
	else
		throw std::runtime_error("Invalid output format: " + f[0] + "\nAllowed values: 0,5,xml,6,tab,100,daa,101,sam,102,103,paf");
}

void init_output()
{
	output_format = unique_ptr<Output_format>(get_output_format());
	if(config.command == Config::view && (output_format->needs_taxon_id_lists || output_format->needs_taxon_nodes || output_format->needs_taxon_scientific_names))
		throw runtime_error("Taxonomy features are not supported for the DAA format.");

	if (*output_format == Output_format::daa && config.multiprocessing)
		throw std::runtime_error("The DAA format is not supported in multiprocessing mode.");
	if (*output_format == Output_format::daa && config.global_ranking_targets)
		throw std::runtime_error("The DAA format is not supported in global ranking mode.");
	if (*output_format == Output_format::taxon && config.toppercent == 100.0 && config.min_bit_score == 0.0)
		config.toppercent = 10.0;
	if (config.toppercent == 100.0) {
		message_stream << "#Target sequences to report alignments for: ";
		if (config.max_alignments == 0) {
			config.max_alignments = std::numeric_limits<size_t>::max();
			message_stream << "unlimited" << endl;
		}
		else
			message_stream << config.max_alignments << endl;
	}
	else
		message_stream << "Percentage range of top alignment score to report hits: " << config.toppercent << endl;
	log_stream << "Format options: transcript=" << bool(output_format->hsp_values & Output::TRANSCRIPT) << " stats=" << bool(output_format->hsp_values & Output::STATS_OR_COORDS) << endl;
}

void Bin1_format::print_query_intro(size_t query_num, const char *query_name, unsigned query_len, TextBuffer &out, bool unaligned, const Search::Config& cfg) const {
	out.write(std::numeric_limits<uint32_t>::max());
	out.write((uint32_t)query_num);
}

void Bin1_format::print_match(const HspContext& r, const Search::Config &metadata, TextBuffer &out) {
	if (r.query_id < r.subject_oid) {
		out.write((uint32_t)r.subject_oid);
		out.write(r.bit_score() / std::max((unsigned)r.query.source().length(), r.subject_len));
	}
}

void Binary_format::print_match(const HspContext& r, const Search::Config& metadata, TextBuffer& out)
{
	out.write((uint32_t)metadata.query->block_id2oid(r.query_id));
	out.write((uint32_t)r.subject_oid);
}


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

#include <limits.h>
#include <algorithm>
#include "data/taxonomy.h"
#include "output_format.h"
#include "util/escape_sequences.h"
#include "output.h"
#include "util/util.h"
#include "run/config.h"
#include "util/sequence/sequence.h"
#include "util/log_stream.h"
#include "util/string/tokenizer.h"

using std::endl;
using std::runtime_error;
using std::string;
using std::vector;

const uint32_t IntermediateRecord::FINISHED = UINT32_MAX;

void IntermediateRecord::read(BinaryBuffer::Iterator& f, const OutputFormat* output_format)
{
	f.read(target_dict_id);
	/*if (config.global_ranking_targets > 0) {
		uint16_t s;
		f.read(s);
		score = s;
		return;
	}*/
	if (*output_format == OutputFormat::daa)
		f.read(target_oid);
	f.read(flag);
	f.read_packed(flag & 3, score);
	f.read(evalue);

	if (output_format->hsp_values == HspValues::NONE)
		return;

	f.read_packed((flag >> 2) & 3, query_begin);
	f.read_varint(query_end);
	f.read_packed((flag >> 4) & 3, subject_begin);

	if (flag_any(output_format->hsp_values, HspValues::TRANSCRIPT))
		transcript.read(f);
	else {
		f.read_varint(subject_end);
		f.read_varint(identities);
		f.read_varint(mismatches);
		f.read_varint(positives);
		f.read_varint(length);
		f.read_varint(gap_openings);
		f.read_varint(gaps);
	}
}

unsigned IntermediateRecord::frame(Loc query_source_len, int align_mode) const {
	if (align_mode == AlignMode::blastx)
		return (flag&(1 << 6)) == 0 ? query_begin % 3 : 3 + (query_source_len - 1 - query_begin) % 3;
	else
		return 0;
}

Interval IntermediateRecord::absolute_query_range() const
{
	if (query_begin < query_end)
		return Interval(query_begin, query_end + 1);
	else
		return Interval(query_end, query_begin + 1);
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

void IntermediateRecord::write(TextBuffer& buf, const Hsp& match, unsigned query_id, DictId target, OId target_oid, const OutputFormat* output_format)
{
	const Interval oriented_range(match.oriented_range());
	buf.write(target);
	if (*output_format == OutputFormat::daa)
		buf.write(target_oid);
	buf.write(get_segment_flag(match));
	buf.write_packed(match.score);
	buf.write(match.evalue);
	if (output_format->hsp_values == HspValues::NONE)
		return;

	buf.write_packed(oriented_range.begin_);
	buf.write_varint(oriented_range.end_);
	buf.write_packed(match.subject_range.begin_);

	if (flag_any(output_format->hsp_values, HspValues::TRANSCRIPT))
		buf << match.transcript.data();
	else {
		buf.write_varint(match.subject_range.end_);
		buf.write_varint(match.identities);
		buf.write_varint(match.mismatches);
		buf.write_varint(match.positives);
		buf.write_varint(match.length);
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

void OutputFormat::print_title(TextBuffer &buf, const char *id, bool full_titles, bool all_titles, const char *separator, const EscapeSequences *esc, bool json_array)
{
	Util::String::Tokenizer<Util::String::StringDelimiters> tok(id, Util::String::StringDelimiters(Util::Seq::FASTA_HEADER_SEP, 2));
	string s;
	int n = 0;
	do {
		if (n++ > 0)
			buf << separator;
		if (json_array)
			buf << '\"';
		tok >> s;
		if (full_titles)
			print_escaped(buf, s, esc);
		else
			print_escaped_until(buf, s.c_str(), Util::Seq::id_delimiters, esc);
		if (json_array)
			buf << '\"';
	} while (tok.good() && all_titles);
}

void print_hsp(Hsp &hsp, const TranslatedSequence &query)
{
	TextBuffer buf;
	//Pairwise_format().print_match(Hsp_context(hsp, 0, query, "", 0, 0, 0, 0, Sequence()), Search::Config(true), buf);
	buf << '\0';
	//std::cout << buf.data() << endl;
}

OutputFormat* get_output_format()
{
	const vector<string> &f = config.output_format;
	if (f.size() == 0) {
		if (config.daa_file == "" || config.command == Config::view)
			return new TabularFormat();
		else if ((config.command == Config::blastp || config.command == Config::blastx) && config.daa_file.length() > 0)
			return new DAA_format();
	}
	if (f[0] == "tab" || f[0] == "6")
		return new TabularFormat();
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
		return new TaxonFormat;
	else if (f[0] == "paf" || f[0] == "103")
		return new PAF_format;
#ifdef WITH_MCL
	else if (f[0] == "bin1")
		return new Bin1_format;
	else if (f[0] == "clus")
		return new Clustering_format(&f[1]);
#endif
	else if (f[0] == "edge")
		return new Output::Format::Edge;
    else if(f[0] == "json-flat" || f[0] == "104")
        return new TabularFormat(true);
	else
		throw std::runtime_error("Invalid output format: " + f[0] + "\nAllowed values: 0,5,xml,6,tab,100,daa,101,sam,102,103,paf");
}

OutputFormat* init_output(int64_t& max_target_seqs)
{
	OutputFormat* output_format = get_output_format();
	if(config.command == Config::view && (output_format->needs_taxon_id_lists || output_format->needs_taxon_nodes || output_format->needs_taxon_scientific_names))
		throw runtime_error("Taxonomy features are not supported for the DAA format.");

	if (*output_format == OutputFormat::daa && config.multiprocessing)
		throw std::runtime_error("The DAA format is not supported in multiprocessing mode.");
	if (*output_format == OutputFormat::daa && config.global_ranking_targets)
		throw std::runtime_error("The DAA format is not supported in global ranking mode.");
	if (*output_format == OutputFormat::taxon && config.toppercent.blank() && config.min_bit_score == 0.0)
		config.toppercent = 10.0;
	if (config.toppercent.present()) {
		if (config.toppercent < 0 || config.toppercent > 100)
			throw runtime_error("Allowed value range for --top is between 0.0 and 100.0");
		if (config.toppercent == 100.0) {
			config.toppercent.unset();
			max_target_seqs = INT64_MAX;
		}
	}
	else {
		max_target_seqs = config.max_target_seqs_.get(DEFAULT_MAX_TARGET_SEQS);
		if (config.max_target_seqs_.present() && config.max_target_seqs_ == 0)
			max_target_seqs = INT64_MAX;
	}
	if (config.toppercent.blank()) {
		message_stream << "#Target sequences to report alignments for: ";
		if (max_target_seqs == INT64_MAX)
			message_stream << "unlimited";
		else
			message_stream << max_target_seqs;
		message_stream << endl;
	}
	else
		message_stream << "Percentage range of top alignment score to report hits: " << config.toppercent << endl;
	if (config.frame_shift != 0 && (output_format->hsp_values != HspValues::NONE || config.query_range_culling))
		output_format->hsp_values = HspValues::TRANSCRIPT;
	log_stream << "DP fields: " << (unsigned)output_format->hsp_values << endl;
	return output_format;
}

namespace Output { namespace Format {

void Edge::print_match(const HspContext& r, Output::Info& info)
{
	info.out.write(Data{ r.query_oid, r.subject_oid, (float)r.qcovhsp(), (float)r.scovhsp(),  r.corrected_bit_score() });
}

}}
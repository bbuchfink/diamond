/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

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

#include "../basic/config.h"
#include "daa_record.h"
#include "../output.h"

DAA_format::DAA_format() :
	OutputFormat(daa, HspValues::TRANSCRIPT, config.salltitles ? Output::Flags::FULL_TITLES : (config.sallseqid ? Output::Flags::ALL_SEQIDS : Output::Flags::NONE))
{}

BinaryBuffer::Iterator DAA_query_record::init(const BinaryBuffer &buf)
{
	BinaryBuffer::Iterator it(buf.begin());
	uint32_t query_len;
	it >> query_len;
	it >> query_name;
	uint8_t flags;
	it >> flags;
	if (file_.mode() == AlignMode::blastp) {
		Packed_sequence seq(it, query_len, false, 5);
		seq.unpack(context[0], 5, query_len);
		query_seq = TranslatedSequence(Sequence(context[0]));
	}
	else {
		const bool have_n = (flags & 1) == 1;
		Packed_sequence seq(it, query_len, have_n, have_n ? 3 : 2);
		seq.unpack(source_seq, have_n ? 3 : 2, query_len);
		translate_query(source_seq, context);
		query_seq = TranslatedSequence(Sequence(source_seq), context);
	}
	return it;
}

void DAA_query_record::Match::read(BinaryBuffer::Iterator &it)
{
	const uint32_t old_subject = subject_id;
	it >> subject_id;
	if (subject_id == old_subject)
		++hsp_num;
	else {
		hsp_num = 0;
		++hit_num;
	}
	uint8_t flag;
	it >> flag;
	it.read_packed(flag & 3, score);
	uint32_t query_begin, subject_begin;
	it.read_packed((flag >> 2) & 3, query_begin);
	it.read_packed((flag >> 4) & 3, subject_begin);
	subject_range.begin_ = (int)subject_begin;
	transcript.read(it);
	subject_name = parent_.file_.ref_name(subject_id);
	subject_len = parent_.file_.ref_len(subject_id);
	if (parent_.file_.mode() == AlignMode::blastx) {
		frame = (flag&(1 << 6)) == 0 ? query_begin % 3 : 3 + (parent_.source_seq.size() - 1 - query_begin) % 3;
		set_translated_query_begin(query_begin, (unsigned)parent_.source_seq.size());
	}
	else if (parent_.file_.mode() == AlignMode::blastp) {
		frame = 0;
		query_range.begin_ = query_begin;
	}
	
	context().parse(nullptr);
	evalue = score_matrix.evalue(score, (Loc)parent_.context[0].size(), subject_len);
	bit_score = score_matrix.bitscore(score);
}

void copy_match_record_raw(BinaryBuffer::Iterator& it, TextBuffer& buf, const std::unordered_map<uint32_t, uint32_t>& subject_map) {
	uint32_t subject_id, query_begin, subject_begin;
	int score;
	it >> subject_id;
	uint8_t flag;
	it >> flag;
	it.read_packed(flag & 3, score);
	it.read_packed((flag >> 2) & 3, query_begin);
	it.read_packed((flag >> 4) & 3, subject_begin);
	buf.write(subject_map.at(subject_id));
	buf.write(flag);
	buf.write_packed(score);
	buf.write_packed(query_begin);
	buf.write_packed(subject_begin);
	const uint8_t t = Packed_operation::terminator();
	uint8_t c;
	do {
		it >> c;
		buf << c;
	} while (c != t);
}

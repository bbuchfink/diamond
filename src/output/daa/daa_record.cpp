/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

#include "basic/config.h"
#include "daa_record.h"
#include "basic/packed_sequence.h"
#include "../output_format.h"
#include "stats/score_matrix.h"

DAAFormat::DAAFormat() :
	OutputFormat(daa, HspValues::TRANSCRIPT, Output::Flags::SSEQID | (config.salltitles ? Output::Flags::FULL_TITLES : (config.sallseqid ? Output::Flags::ALL_SEQIDS : Output::Flags::NONE)))
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
		PackedSequence seq(it, query_len, false, 5);
		seq.unpack(context[0], 5, query_len);
		query_seq = TranslatedSequence(Sequence(context[0]));
	}
	else {
		const bool have_n = (flags & 1) == 1;
		PackedSequence seq(it, query_len, have_n, have_n ? 3 : 2);
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
	
	*(Hsp*)this = context().parse(nullptr).hsp();
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
	const uint8_t t = PackedOperation::terminator();
	uint8_t c;
	do {
		it >> c;
		buf << c;
	} while (c != t);
}

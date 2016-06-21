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

#include "daa_record.h"

const string & DAA_query_record::Match::query_name() const
{
	return parent_.query_name;
}

const vector<Letter>& DAA_query_record::Match::query() const
{
	return parent_.context[frame];
}

uint64_t DAA_query_record::Match::db_letters() const
{
	return parent_.file_.db_letters();
}

unsigned DAA_query_record::Match::query_end() const
{
	if (parent_.file_.mode() == mode_blastp) {
		return query_begin + translated_query_len - 1;
	}
	else if (parent_.file_.mode() == mode_blastx) {
		int len = (int)translated_query_len * 3 * (frame>2 ? -1 : 1);
		return (int)query_begin + (len > 0 ? -1 : 1) + len;
	}
	else if (parent_.file_.mode() == mode_blastn) {
		int len = (int)translated_query_len*(frame>0 ? -1 : 1);
		return (int)query_begin + (len > 0 ? -1 : 1) + len;
	}
	else
		return 0;
}

void DAA_query_record::Match::parse()
{
	len = identities = mismatches = gap_openings = positives = gaps = 0;
	unsigned d = 0;
	Position_iterator i = begin();

	for (; i.good(); ++i) {
		++len;
		switch (i.op()) {
		case op_match:
			++identities;
			++positives;
			d = 0;
			break;
		case op_substitution:
			++mismatches;
			if (i.score() > 0)
				++positives;
			d = 0;
			break;
		case op_insertion:
		case op_deletion:
			if (d == 0)
				++gap_openings;
			++d;
			++gaps;
			break;
		}
	}

	translated_query_len = i.query_pos - translated_query_begin;
	subject_len = i.subject_pos - subject_begin;
}

Binary_buffer::Iterator DAA_query_record::init(const Binary_buffer &buf)
{
	Binary_buffer::Iterator it(buf.begin());
	uint32_t query_len;
	it >> query_len;
	it >> query_name;
	uint8_t flags;
	it >> flags;
	if (file_.mode() == mode_blastp) {
		Packed_sequence seq(it, query_len, false, 5);
		seq.unpack(context[0], 5, query_len);
	}
	else {
		const bool have_n = (flags & 1) == 1;
		Packed_sequence seq(it, query_len, have_n, have_n ? 3 : 2);
		seq.unpack(source_seq, have_n ? 3 : 2, query_len);
		translate_query(source_seq, context);
	}
	return it;
}

Binary_buffer::Iterator& operator>>(Binary_buffer::Iterator &it, DAA_query_record::Match &r)
{
	const uint32_t old_subject = r.subject_id;
	it >> r.subject_id;
	if (r.subject_id == old_subject)
		++r.hsp_num;
	else {
		r.hsp_num = 0;
		++r.hit_num;
	}
	uint8_t flag;
	it >> flag;
	it.read_packed(flag & 3, r.score);
	it.read_packed((flag >> 2) & 3, r.query_begin);
	it.read_packed((flag >> 4) & 3, r.subject_begin);
	r.transcript.read(it);
	r.subject_name = r.parent_.file_.ref_name(r.subject_id);
	r.total_subject_len = r.parent_.file_.ref_len(r.subject_id);
	if (r.parent_.file_.mode() == mode_blastx) {
		r.frame = (flag&(1 << 6)) == 0 ? r.query_begin % 3 : 3 + (r.parent_.source_seq.size() - 1 - r.query_begin) % 3;
		r.translated_query_begin = query_translated_begin(r.query_begin, r.frame, (unsigned)r.parent_.source_seq.size(), true);
	}
	else if (r.parent_.file_.mode() == mode_blastp) {
		r.frame = 0;
		r.translated_query_begin = r.query_begin;
	}
	else {
		r.frame = (flag&(1 << 6)) == 0 ? 0 : 1;
		r.translated_query_begin = query_translated_begin(r.query_begin, r.frame, (unsigned)r.parent_.source_seq.size(), false);
	}
	r.parse();
	return it;
}
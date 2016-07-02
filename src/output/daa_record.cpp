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
#include "output.h"

Binary_buffer::Iterator DAA_query_record::init(const Binary_buffer &buf)
{
	Binary_buffer::Iterator it(buf.begin());
	uint32_t query_len;
	it >> query_len;
	it >> query_name;
	uint8_t flags;
	it >> flags;
	if (file_.mode() == Align_mode::blastp) {
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
	uint32_t query_begin;
	it.read_packed((flag >> 2) & 3, query_begin);
	it.read_packed((flag >> 4) & 3, r.subject_range.begin_);
	r.transcript.read(it);
	r.subject_name = r.parent_.file_.ref_name(r.subject_id);
	r.subject_len = r.parent_.file_.ref_len(r.subject_id);
	if (r.parent_.file_.mode() == Align_mode::blastx) {
		r.frame = (flag&(1 << 6)) == 0 ? query_begin % 3 : 3 + (r.parent_.source_seq.size() - 1 - query_begin) % 3;
		r.set_translated_query_begin(query_begin, (unsigned)r.parent_.source_seq.size());
	}
	else if (r.parent_.file_.mode() == Align_mode::blastp) {
		r.frame = 0;
		r.query_range.begin_ = query_begin;
	}
	r.context().parse().set_query_source_range(query_begin);
	return it;
}

Hsp_data::Hsp_data(const Intermediate_record &r, unsigned query_source_len):
	score(r.score),
	transcript(r.transcript)
{
	subject_range.begin_ = r.subject_begin;
	if (align_mode.mode == Align_mode::blastx) {
		frame = (r.flag&(1 << 6)) == 0 ? r.query_begin % 3 : 3 + (query_source_len - 1 - r.query_begin) % 3;
		set_translated_query_begin(r.query_begin, query_source_len);
	}
	else {
		frame = 0;
		query_range.begin_ = r.query_begin;
	}
}
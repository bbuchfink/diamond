/****
Copyright (c) 2014, University of Tuebingen
Author: Benjamin Buchfink
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

#ifndef DAA_WRITE_H_
#define DAA_WRITE_H_

#include "daa_file.h"

struct Intermediate_record
{
	void read(Buffered_file &f)
	{
		f.read(query_id);
		f.read(subject_id);
		f.read(flag);
		f.read_packed(flag & 3, score);
		f.read_packed((flag>>2)&3, query_begin);
		f.read_packed((flag>>4)&3, subject_begin);
		transcript.read(f);
	}
	uint32_t query_id, subject_id, score, query_begin, subject_begin;
	uint8_t flag;
	Packed_transcript transcript;
};

struct DAA_output
{

	DAA_output():
		f_ (program_options::daa_file),
		h2_ (ref_header.sequences,
				ref_header.letters,
				program_options::gap_open,
				program_options::gap_extend,
				program_options::reward,
				program_options::penalty,
				score_matrix::get().k(),
				score_matrix::get().lambda(),
				program_options::matrix,
				(Align_mode)program_options::command)
	{
		DAA_header1 h1;
		f_.write(&h1, 1);
		h2_.block_type[0] = DAA_header2::alignments;
		h2_.block_type[1] = DAA_header2::ref_names;
		h2_.block_type[2] = DAA_header2::ref_lengths;
		f_.write(&h2_, 1);
	}

	template<typename _val>
	static void write_query_record(Text_buffer &buf, const sequence<const char> &query_name, const sequence<const _val> &query)
	{
		buf.write((uint32_t)0);
		uint32_t l = query.length();
		buf.write(l);
		buf.write_c_str(query_name.c_str(), find_first_of(query_name.c_str(), Const::id_delimiters));
		Packed_sequence s (query);
		uint8_t flags = s.has_n() ? 1 : 0;
		buf.write(flags);
		buf << s.data();
	}

	static void write_record(Text_buffer &buf, const Intermediate_record &r)
	{
		buf.write(r.subject_id).write(r.flag);
		buf.write_packed(r.score);
		buf.write_packed(r.query_begin);
		buf.write_packed(r.subject_begin);
		buf << r.transcript.data();
	}

	template<typename _val>
	static void write_record(Text_buffer &buf,
			const Segment<_val> &match,
			size_t query_source_len,
			const sequence<const _val> &query,
			unsigned query_id,
			const vector<char> &transcript_buf)
	{
		buf.write(ref_map.get<_val>(current_ref_block, match.subject_id_));
		buf.write(get_segment_flag(match));
		buf.write_packed(match.score_);
		buf.write_packed(match.traceback_->query_begin_);
		buf.write_packed(match.traceback_->subject_begin_);
		const unsigned qbegin = query_translated_begin<_val>(match.traceback_->query_begin_, match.frame_, query_source_len, query_translated());
		print_packed(match.traceback_->transcript_right_, match.traceback_->transcript_left_, transcript_buf, buf, query, ref_seqs<_val>::get()[match.subject_id_], qbegin, match.traceback_->subject_begin_);
	}

	void finish()
	{
		uint32_t size = 0;
		f_.write(&size, 1);
		h2_.block_size[0] = f_.tell() - sizeof(DAA_header1) - sizeof(DAA_header2);
		h2_.db_seqs_used = ref_map.next_;
		h2_.query_records = statistics.get(Statistics::ALIGNED);

		size_t s = 0;
		for(ptr_vector<string>::const_iterator i = ref_map.name_.begin(); i != ref_map.name_.end(); ++i) {
			f_.write_c_str(*i);
			s += i->length()+1;
		}
		h2_.block_size[1] = s;

		f_.write(ref_map.len_, false);
		h2_.block_size[2] = ref_map.len_.size() * sizeof(uint32_t);

		f_.seek(sizeof(DAA_header1));
		f_.write(&h2_, 1);

		f_.close();
	}

	Output_stream& stream()
	{ return f_; }

private:

	Output_stream f_;
	DAA_header2 h2_;

};

#endif /* DAA_WRITE_H_ */

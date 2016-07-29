/****
Copyright (c) 2014-2016, University of Tuebingen, Benjamin Buchfink
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

#include <limits>
#include <stdint.h>
#include "output.h"
#include "daa_file.h"
#include "../basic/score_matrix.h"
#include "../data/reference.h"
#include "../align/align.h"
#include "../basic/packed_sequence.h"

inline void init_daa(Output_stream &f)
{
	DAA_header1 h1;
	f.typed_write(&h1, 1);	
	DAA_header2 h2_;
	f.typed_write(&h2_, 1);
}

inline size_t write_daa_query_record(Text_buffer &buf, const char *query_name, const sequence &query)
{
	size_t seek_pos = buf.size();
	buf.write((uint32_t)0);
	uint32_t l = (uint32_t)query.length();
	buf.write(l);
	buf.write_c_str(query_name, find_first_of(query_name, Const::id_delimiters));
	Packed_sequence s(query, align_mode.input_sequence_type);
	uint8_t flags = s.has_n() ? 1 : 0;
	buf.write(flags);
	buf << s.data();
	return seek_pos;
}

inline void finish_daa_query_record(Text_buffer &buf, size_t seek_pos)
{
	*(uint32_t*)(&buf[seek_pos]) = (uint32_t)(buf.size() - seek_pos - sizeof(uint32_t));
}

inline void write_daa_record(Text_buffer &buf, const Intermediate_record &r)
{
	buf.write(r.subject_id).write(r.flag);
	buf.write_packed(r.score);
	buf.write_packed(r.query_begin);
	buf.write_packed(r.subject_begin);
	buf << r.transcript.data();
}

inline void write_daa_record(Text_buffer &buf, const Hsp_data &match, unsigned query_id, unsigned subject_id)
{
	buf.write(ref_map.get(current_ref_block, subject_id));
	buf.write(get_segment_flag(match));
	buf.write_packed(match.score);
	buf.write_packed(match.oriented_range().begin_);
	buf.write_packed(match.subject_range.begin_);
	buf << match.transcript.data();
}

inline void finish_daa(Output_stream &f)
{
	DAA_header2 h2_(ref_header.sequences,
		config.db_size,
		config.gap_open,
		config.gap_extend,
		config.reward,
		config.penalty,
		score_matrix.k(),
		score_matrix.lambda(),
		config.max_evalue,
		to_lower_case(config.matrix),
		align_mode.mode);

	h2_.block_type[0] = DAA_header2::alignments;
	h2_.block_type[1] = DAA_header2::ref_names;
	h2_.block_type[2] = DAA_header2::ref_lengths;

	uint32_t size = 0;
	f.typed_write(&size, 1);
	h2_.block_size[0] = f.tell() - sizeof(DAA_header1) - sizeof(DAA_header2);
	h2_.db_seqs_used = ref_map.next_;
	h2_.query_records = statistics.get(Statistics::ALIGNED);

	size_t s = 0;
	for (Ptr_vector<string>::const_iterator i = ref_map.name_.begin(); i != ref_map.name_.end(); ++i) {
		f.write_c_str((*i)->c_str());
		s += (*i)->length() + 1;
	}
	h2_.block_size[1] = s;

	f.write(ref_map.len_, false);
	h2_.block_size[2] = ref_map.len_.size() * sizeof(uint32_t);

	f.seekp(sizeof(DAA_header1));
	f.typed_write(&h2_, 1);
}

#endif /* DAA_WRITE_H_ */

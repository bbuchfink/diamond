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

#ifndef OUTPUT_H_
#define OUTPUT_H_

#include "../util/binary_file.h"
#include "../basic/packed_transcript.h"
#include "../util/text_buffer.h"
#include "../align/align_struct.h"
#include "../data/reference.h"

inline unsigned get_length_flag(unsigned x)
{
	if (x <= (unsigned)std::numeric_limits<uint8_t>::max())
		return 0;
	else if (x <= (unsigned)std::numeric_limits<uint16_t>::max())
		return 1;
	return 2;
}

inline unsigned get_rev_flag(unsigned frame)
{
	return frame > 2 ? 1 : 0;
}

inline uint8_t get_segment_flag(const Hsp_data &match)
{
	unsigned rev = get_rev_flag(match.frame);
	return (uint8_t)(get_length_flag(match.score)
		| (get_length_flag(match.oriented_range().begin_) << 2)
		| (get_length_flag(match.subject_range.begin_) << 4)
		| rev << 6);
}

inline uint8_t get_segment_flag(const Hsp_context &match)
{
	unsigned rev = get_rev_flag(match.frame());
	return (uint8_t)(get_length_flag(match.score())
		| (get_length_flag(match.oriented_query_range().begin_) << 2)
		| (get_length_flag(match.subject_range().begin_) << 4)
		| rev << 6);
}

struct Intermediate_record
{
	void read(Buffered_file &f)
	{
		f.read(query_id);
		f.read(subject_id);
		f.read(flag);
		f.read_packed(flag & 3, score);
		f.read_packed((flag >> 2) & 3, query_begin);
		f.read_packed((flag >> 4) & 3, subject_begin);
		transcript.read(f);
	}
	void read(Binary_buffer::Iterator &f)
	{
		f.read(subject_id);
		f.read(flag);
		f.read_packed(flag & 3, score);
		f.read_packed((flag >> 2) & 3, query_begin);
		f.read_packed((flag >> 4) & 3, subject_begin);
		transcript.read(f);
	}
	static size_t write_query_intro(Text_buffer &buf, unsigned query_id)
	{
		size_t seek_pos = buf.size();
		buf.write(query_id).write(0u);
		return seek_pos;
	}
	static void finish_query(Text_buffer &buf, size_t seek_pos)
	{
		*(unsigned*)(&buf[seek_pos + sizeof(unsigned)]) = safe_cast<unsigned>(buf.size() - seek_pos - sizeof(unsigned) * 2);
	}
	static void write(Text_buffer &buf, const Hsp_data &match, unsigned query_id, unsigned subject_id)
	{
		buf.write(ref_map.get(current_ref_block, subject_id))
			.write(get_segment_flag(match))
			.write_packed(match.score)
			.write_packed(match.oriented_range().begin_)
			.write_packed(match.subject_range.begin_)
			<< match.transcript.data();
	}
	static void finish_file(Output_stream &f)
	{
		unsigned x = finished;
		f.typed_write(&x, 1);
	}
	enum { finished = 0xffffffffu };
	uint32_t query_id, subject_id, score, query_begin, subject_begin;
	uint8_t flag;
	Packed_transcript transcript;
};

#ifndef ST_JOIN
void join_blocks(unsigned ref_blocks, Output_stream &master_out, const vector<Temp_file> &tmp_file);
#endif

#endif
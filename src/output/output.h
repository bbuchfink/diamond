/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#ifndef OUTPUT_H_
#define OUTPUT_H_

#include <map>
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
		f.write(&x, 1);
	}
	enum { finished = 0xffffffffu };
	uint32_t query_id, subject_id, score, query_begin, subject_begin;
	uint8_t flag;
	Packed_transcript transcript;
};

void join_blocks(unsigned ref_blocks, Output_stream &master_out, const vector<Temp_file> &tmp_file);

struct Output_sink
{
	Output_sink(size_t begin, Output_stream *f) :
		f_(f),
		next_(begin),
		size_(0),
		max_size_(0)
	{}
	void push(size_t n, Text_buffer *buf);
	size_t size() const
	{
		return size_;
	}
	size_t max_size() const
	{
		return max_size_;
	}
	static Output_sink& get()
	{
		return *instance;
	}
	size_t next() const
	{
		return next_;
	}
	static auto_ptr<Output_sink> instance;
private:
	void flush(Text_buffer *buf);
	tthread::mutex mtx_;
	Output_stream* const f_;
	std::map<size_t, Text_buffer*> backlog_;
	size_t next_, size_, max_size_;
};

void heartbeat_worker(size_t qend);

#endif
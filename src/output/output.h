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

#ifndef OUTPUT_H_
#define OUTPUT_H_

#include <memory>
#include <map>
#include <mutex>
#include "../util/io/output_file.h"
#include "../basic/packed_transcript.h"
#include "../util/text_buffer.h"
#include "../data/reference.h"
#include "../basic/match.h"
#include "../data/ref_dictionary.h"
#include "../basic/parameters.h"
#include "../data/metadata.h"
#include "../util/io/consumer.h"

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

inline uint8_t get_segment_flag(const Hsp &match)
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

struct IntermediateRecord
{
	void read(BinaryBuffer::Iterator &f)
	{
		f.read(subject_id);
		f.read(flag);
		f.read_packed(flag & 3, score);
		f.read_packed((flag >> 2) & 3, query_begin);
		f.read_varint(query_end);
		f.read_packed((flag >> 4) & 3, subject_begin);
		transcript.read(f);
	}
	interval absolute_query_range() const
	{
		if (query_begin < query_end)
			return interval(query_begin, query_end + 1);
		else
			return interval(query_end, query_begin + 1);
	}
	static size_t write_query_intro(TextBuffer &buf, unsigned query_id)
	{
		size_t seek_pos = buf.size();
		buf.write(query_id).write(0u);
		return seek_pos;
	}
	static void finish_query(TextBuffer &buf, size_t seek_pos)
	{
		*(unsigned*)(&buf[seek_pos + sizeof(unsigned)]) = safe_cast<unsigned>(buf.size() - seek_pos - sizeof(unsigned) * 2);
	}
	static void write(TextBuffer &buf, const Hsp &match, unsigned query_id, size_t subject_id)
	{
		const interval oriented_range (match.oriented_range());
		buf.write(ReferenceDictionary::get().get(current_ref_block, subject_id))
			.write(get_segment_flag(match))
			.write_packed(match.score)
			.write_packed(oriented_range.begin_)
			.write_varint(oriented_range.end_)
			.write_packed(match.subject_range.begin_)
			<< match.transcript.data();
	}
	static void finish_file(Consumer &f)
	{
		const uint32_t i = FINISHED;
		f.consume(reinterpret_cast<const char*>(&i), 4);
	}
	static const uint32_t FINISHED = 0xffffffffu;
	uint32_t query_id, subject_id, score, query_begin, subject_begin, query_end;
	uint8_t flag;
	Packed_transcript transcript;
};

void join_blocks(unsigned ref_blocks, Consumer &master_out, const PtrVector<TempFile> &tmp_file, const Parameters &params, const Metadata &metadata, DatabaseFile &db_file);

struct OutputSink
{
	OutputSink(size_t begin, Consumer *f) :
		f_(f),
		begin_(begin),
		next_(begin),
		size_(0),
		max_size_(0)
	{}
	void push(size_t n, TextBuffer *buf);
	size_t size() const
	{
		return size_;
	}
	size_t max_size() const
	{
		return max_size_;
	}
	static OutputSink& get()
	{
		return *instance;
	}
	size_t next() const
	{
		return next_;
	}
	size_t begin() const {
		return begin_;
	}
	static std::unique_ptr<OutputSink> instance;
private:
	void flush(TextBuffer *buf);
	std::mutex mtx_;
	Consumer* const f_;
	std::map<size_t, TextBuffer*> backlog_;
	size_t begin_, next_, size_, max_size_;
};

void heartbeat_worker(size_t qend);

#endif
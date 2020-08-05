/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
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

#pragma once
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
#include "output_format.h"

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

		if (!output_format->needs_transcript && !output_format->needs_stats)
			return;

		f.read_packed((flag >> 2) & 3, query_begin);
		f.read_varint(query_end);
		f.read_packed((flag >> 4) & 3, subject_begin);

		if (output_format->needs_transcript)
			transcript.read(f);
		else if (output_format->needs_stats) {
			f.read_varint(subject_end);
			f.read_varint(identities);
			f.read_varint(mismatches);
			f.read_varint(positives);
			f.read_varint(gap_openings);
			f.read_varint(gaps);
		}
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
		buf.write((uint32_t)query_id).write((uint32_t)0);
		return seek_pos;
	}
	static void finish_query(TextBuffer &buf, size_t seek_pos)
	{
		*(uint32_t*)(&buf[seek_pos + sizeof(uint32_t)]) = safe_cast<uint32_t>(buf.size() - seek_pos - sizeof(uint32_t) * 2);
	}
	static void write(TextBuffer &buf, const Hsp &match, unsigned query_id, size_t subject_id)
	{
		const interval oriented_range (match.oriented_range());
		buf.write(ReferenceDictionary::get().get(current_ref_block, subject_id));
		buf.write(get_segment_flag(match));
		buf.write_packed(match.score);
		if (!output_format->needs_transcript && !output_format->needs_stats)
			return;

		buf.write_packed(oriented_range.begin_);
		buf.write_varint(oriented_range.end_);
		buf.write_packed(match.subject_range.begin_);

		if(output_format->needs_transcript)
			buf << match.transcript.data();
		else if (output_format->needs_stats) {
			buf.write_varint(match.subject_range.end_);
			buf.write_varint(match.identities);
			buf.write_varint(match.mismatches);
			buf.write_varint(match.positives); 
			buf.write_varint(match.gap_openings);
			buf.write_varint(match.gaps);
		}
	}
	static void finish_file(Consumer &f)
	{
		const uint32_t i = FINISHED;
		f.consume(reinterpret_cast<const char*>(&i), 4);
	}
	static const uint32_t FINISHED = 0xffffffffu;
	uint32_t query_id, subject_id, score, query_begin, subject_begin, query_end, subject_end, identities, mismatches, positives, gap_openings, gaps;
	uint8_t flag;
	Packed_transcript transcript;
};

void join_blocks(unsigned ref_blocks, Consumer &master_out, const PtrVector<TempFile> &tmp_file, const Parameters &params, const Metadata &metadata, DatabaseFile &db_file,
					const vector<string> tmp_file_names = vector<string>());

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

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
#include <stdint.h>
#include "../util/io/output_file.h"
#include "../basic/packed_transcript.h"
#include "../util/text_buffer.h"
#include "../data/reference.h"
#include "../basic/match.h"
#include "../util/io/consumer.h"
#include "output_format.h"
#include "../run/config.h"
#include "../util/data_structures/reorder_queue.h"

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

inline uint8_t get_segment_flag(const HspContext &match)
{
	unsigned rev = get_rev_flag(match.frame());
	return (uint8_t)(get_length_flag(match.score())
		| (get_length_flag(match.oriented_query_range().begin_) << 2)
		| (get_length_flag(match.subject_range().begin_) << 4)
		| rev << 6);
}

struct IntermediateRecord
{
	void read(BinaryBuffer::Iterator& f, const OutputFormat* output_format);
	unsigned frame(Loc query_source_len, int align_mode) const;
	Interval absolute_query_range() const;
	static size_t write_query_intro(TextBuffer& buf, unsigned query_id);
	static void finish_query(TextBuffer& buf, size_t seek_pos);
	static void write(TextBuffer& buf, const Hsp& match, unsigned query_id, DictId target, OId target_oid, const OutputFormat* output_format);
	static void write(TextBuffer& buf, uint32_t target_block_id, int score, const Search::Config& cfg);
	static void finish_file(Consumer& f);
	static bool stats_mode(const HspValues v) {
		return !flag_any(v, HspValues::TRANSCRIPT) && v != HspValues::NONE;
	}
	static const uint32_t FINISHED;
	BlockId query_id;
	DictId target_dict_id;
	OId target_oid;
	uint32_t score, query_begin, subject_begin, query_end, subject_end, identities, mismatches, positives, length, gap_openings, gaps;
	double evalue;
	uint8_t flag;
	Packed_transcript transcript;
};

void join_blocks(int64_t ref_blocks, Consumer &master_out, const PtrVector<TempFile> &tmp_file, Search::Config& cfg, SequenceFile &db_file,
					const std::vector<std::string> tmp_file_names = std::vector<std::string>());

struct OutputWriter {
    OutputWriter(Consumer* file_, char sep = char(0), bool first = true):
    file_(file_),
    sep(sep),
    first(first)
    {};
	void operator()(TextBuffer* buf) {
        if(!first && sep != char(0))
        {
            file_ ->consume(&sep, 1);
        }
		file_->consume(buf->data(), buf->size());
        first = false;
	}
	Consumer* file_;
    bool first;
    char sep;
};

extern std::unique_ptr<ReorderQueue<TextBuffer*, OutputWriter>> output_sink;

void heartbeat_worker(size_t qend, const Search::Config* cfg);

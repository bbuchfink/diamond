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

#pragma once
#include <memory>
#include <stdint.h>
#include "basic/packed_transcript.h"
#include "util/text_buffer.h"
#include "basic/match.h"
#include "util/io/consumer.h"
#include "output_format.h"
#include "run/config.h"
#include "util/data_structures/reorder_queue.h"
#include "util/ptr_vector.h"

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
	PackedTranscript transcript;
};

void join_blocks(int64_t ref_blocks, Consumer &master_out, const PtrVector<TempFile> &tmp_file, Search::Config& cfg, SequenceFile &db_file,
					const std::vector<std::string> tmp_file_names = std::vector<std::string>());

struct OutputWriter {
	OutputWriter(Consumer* file_, char sep = char(0), bool first = true) :
		file_(file_),
		first(first),
		sep(sep)
	{};
	void operator()(TextBuffer* buf) {
        if(!first && sep != '\0')
        {
            file_->consume(&sep, 1);
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
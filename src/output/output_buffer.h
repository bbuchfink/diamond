/****
Copyright (c) 2014, University of Tuebingen
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****
Author: Benjamin Buchfink
****/

#ifndef OUTPUT_BUFFER_H_
#define OUTPUT_BUFFER_H_

#include "../util/text_buffer.h"
#include "output_format.h"
#include "daa_write.h"

unsigned get_length_flag(unsigned x)
{
	if(x <= (unsigned)std::numeric_limits<uint8_t>::max())
		return 0;
	else if(x <= (unsigned)std::numeric_limits<uint16_t>::max())
		return 1;
	else
		return 2;
}

template<typename _val>
unsigned get_rev_flag(unsigned frame)
{ return frame > 0 ? 1 : 0; }

template<>
unsigned get_rev_flag<Amino_acid>(unsigned frame)
{ return frame > 2 ? 1 : 0; }

template<typename _val>
uint8_t get_segment_flag(const Segment<_val> &match)
{
	unsigned rev = get_rev_flag<_val>(match.frame_);
	return (uint8_t)(get_length_flag(match.score_)
				| (get_length_flag(match.traceback_->query_begin_)<<2)
				| (get_length_flag(match.traceback_->subject_begin_)<<4)
				| rev<<6);
}

template<typename _val>
void write_intermediate_record(Text_buffer &buf,
			const Segment<_val> &match,
			size_t query_source_len,
			const sequence<const _val> &query,
			unsigned query_id)
{
	buf.write(query_id)
		.write(ref_map.get<_val>(current_ref_block, match.subject_id_))
		.write(get_segment_flag(match))
		.write_packed(match.score_)
		.write_packed(match.traceback_->query_begin_)
		.write_packed(match.traceback_->subject_begin_);
	const unsigned qbegin = query_translated_begin<_val>(match.traceback_->query_begin_, match.frame_, query_source_len);
	match.traceback_->transcript_->print_packed(buf, query, ref_seqs<_val>::get()[match.subject_id_], qbegin, match.traceback_->subject_begin_);
}

template<typename _val>
struct Output_buffer : public Text_buffer
{
	virtual void print_match(const Segment<_val> &match,
		size_t query_source_len,
		const sequence<const _val> &query,
		unsigned query_id)
	{ DAA_output::write_record(*this, match, query_source_len, query, query_id); }
	virtual void write_query_record(unsigned query_id)
	{
		query_begin_ = this->size();
		if(query_translated())
			DAA_output::write_query_record(*this, query_ids::get()[query_id], query_source_seqs::get()[query_id]);
		else
			DAA_output::write_query_record(*this, query_ids::get()[query_id], query_seqs<_val>::get()[query_id]);
	}
	virtual void finish_query_record()
	{ *(uint32_t*)(this->data_+query_begin_) = this->size() - query_begin_ - sizeof(uint32_t); }
	virtual ~Output_buffer()
	{ }
private:
	size_t query_begin_;
};

template<typename _val>
struct Temp_output_buffer : public Output_buffer<_val>
{
	virtual void print_match(const Segment<_val> &match,
				size_t query_source_len,
				const sequence<const _val> &query,
				unsigned query_id)
	{ write_intermediate_record(*this, match, query_source_len, query, query_id); }
	virtual void write_query_record(unsigned query_id)
	{ }
	virtual void finish_query_record()
	{ }
	virtual ~Temp_output_buffer()
	{ }
};

#endif /* OUTPUT_BUFFER_H_ */

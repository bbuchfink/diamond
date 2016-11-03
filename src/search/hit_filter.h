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

#ifndef HIT_FILTER_H_
#define HIT_FILTER_H_

#include <vector>
#include <limits>
#include "trace_pt_buffer.h"
#include "../dp/smith_waterman.h"
#include "../basic/sequence.h"
#include "../data/queries.h"
#include "../data/reference.h"

using std::vector;

#ifdef __SSE2__

struct hit_filter
{

	hit_filter(Statistics &stats,
			   Loc q_pos,
			   Trace_pt_buffer::Iterator &out):
		q_num_ (std::numeric_limits<unsigned>::max()),
		seed_offset_ (std::numeric_limits<unsigned>::max()),
		stats_ (stats),
		q_pos_ (q_pos),
		out_ (out),
		subjects_ (TLS::get(subjects_ptr))
	{ subjects_.clear(); }

	void push(Loc subject, int score)
	{
		if(score >= config.min_hit_raw_score)
			push_hit(subject);
		else
			subjects_.push_back(ref_seqs::data_->fixed_window_infix(subject+ config.seed_anchor));
	}

	void finish()
	{
		if(subjects_.size() == 0)
			return;
		unsigned left;
		sequence query (query_seqs::data_->window_infix(q_pos_ + config.seed_anchor, left));
		smith_waterman(query,
				subjects_,
				config.hit_band,
				left,
				config.gap_open + config.gap_extend,
				config.gap_extend,
				config.min_hit_raw_score,
				*this,
				uint8_t(),
				stats_);
	}

	void push_hit(Loc subject)
	{
		if(q_num_ == std::numeric_limits<unsigned>::max()) {
			std::pair<size_t,size_t> l (query_seqs::data_->local_position(q_pos_));
			q_num_ = (unsigned)l.first;
			seed_offset_ = (unsigned)l.second;
		}
		assert(subject < ref_seqs::get().raw_len());
		out_.push(hit  (q_num_, subject, seed_offset_));
		stats_.inc(Statistics::TENTATIVE_MATCHES4);
	}

	void operator()(int i, const sequence &seq, int score)
	{ push_hit(ref_seqs::data_->position(seq.data()+config.window- config.seed_anchor)); stats_.inc(Statistics::GAPPED_HITS); }

private:

	unsigned q_num_, seed_offset_;
	Statistics  &stats_;
	Loc q_pos_;
	Trace_pt_buffer::Iterator &out_;
	vector<sequence> &subjects_;
	
	static TLS_PTR vector<sequence> *subjects_ptr;

};

#endif

#endif /* HIT_FILTER_H_*/

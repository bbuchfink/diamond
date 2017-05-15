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
				score_matrix.gap_open() + score_matrix.gap_extend(),
				score_matrix.gap_extend(),
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

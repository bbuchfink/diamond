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
#include <boost/thread/tss.hpp>
#include "trace_pt_buffer.h"
#include "../dp/smith_waterman.h"
#include "../basic/sequence.h"

using std::vector;
using boost::thread_specific_ptr;

template<typename _val, typename _locr, typename _locq, typename _locl>
struct hit_filter
{

	hit_filter(Statistics &stats,
			   _locq q_pos,
			   typename Trace_pt_buffer<_locr,_locl>::Iterator &out):
		q_num_ (std::numeric_limits<unsigned>::max()),
		seed_offset_ (std::numeric_limits<unsigned>::max()),
		stats_ (stats),
		q_pos_ (q_pos),
		out_ (out),
		subjects_ (&s2)
		//subjects_ (subjects_ptr)
	{ subjects_->clear(); }

	void push(_locr subject, int score)
	{
		if(score >= program_options::min_hit_score)
			push_hit(subject);
		else
			subjects_->push_back(ref_seqs<_val>::data_->fixed_window_infix(subject+Const::seed_anchor));
	}

	void finish()
	{
		if(subjects_->size() == 0)
			return;
		unsigned left;
		sequence<const _val> query (query_seqs<_val>::data_->window_infix(q_pos_ + Const::seed_anchor, left));
		smith_waterman(query,
				*subjects_,
				program_options::hit_band,
				left,
				program_options::gap_open + program_options::gap_extend,
				program_options::gap_extend,
				program_options::min_hit_score,
				*this,
				uint8_t(),
				stats_);
	}

	void push_hit(_locr subject)
	{
		if(q_num_ == std::numeric_limits<unsigned>::max()) {
			std::pair<size_t,size_t> l (query_seqs<_val>::data_->local_position(q_pos_));
			q_num_ = l.first;
			seed_offset_ = l.second;
		}
		//cout << "query=" << q_num_ << " so=" << seed_offset_ << " subject=" << subject << endl;
		assert(subject < ref_seqs<_val>::get().raw_len());
		out_.push(hit<_locr,_locl> (q_num_, subject, seed_offset_));
		stats_.inc(Statistics::TENTATIVE_MATCHES3);
	}

	void operator()(int i, const sequence<const _val> &seq, int score)
	{ push_hit(ref_seqs<_val>::data_->position(seq.data()+program_options::window-Const::seed_anchor)); stats_.inc(Statistics::GAPPED_HITS); }

private:

	unsigned q_num_, seed_offset_;
	Statistics  &stats_;
	_locq q_pos_;
	typename Trace_pt_buffer<_locr,_locl>::Iterator &out_;
	//Tls<vector<sequence<const _val> > > subjects_;
	vector<sequence<const _val> > s2;
	vector<sequence<const _val> >* subjects_;

	//static thread_specific_ptr<vector<sequence<const _val> > > subjects_ptr;

};

/*template<typename _val, typename _locr, typename _locq, typename _locl>
thread_specific_ptr<vector<sequence<const _val> > > hit_filter<_val,_locr,_locq,_locl>::subjects_ptr;*/

#endif /* HIT_FILTER_H_*/

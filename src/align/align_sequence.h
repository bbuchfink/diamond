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

#ifndef ALIGN_SEQUENCE_H_
#define ALIGN_SEQUENCE_H_

#include <vector>
#include "../dp/floating_sw.h"

using std::vector;

template<typename _val, typename _locr, typename _locl>
void align_sequence(vector<Segment<_val> > &matches,
		Statistics &stat,
		vector<local_match<_val> > &local,
		unsigned *padding,
		size_t db_letters,
		unsigned dna_len,
		typename Trace_pt_buffer<_locr,_locl>::Vector::iterator &begin,
		typename Trace_pt_buffer<_locr,_locl>::Vector::iterator &end,
		vector<char> &transcript_buf)
{
	std::sort(begin, end, hit<_locr,_locl>::cmp_normalized_subject);
	const unsigned q_num (begin->query_);
	const sequence<const _val> query (query_seqs<_val>::get()[q_num]);
	const unsigned frame = q_num % query_contexts();
	const unsigned query_len = query.length();
	padding[frame] = program_options::read_padding<_val>(query_len);

	const Sequence_set<_val> *ref = ref_seqs<_val>::data_;
	for(typename Trace_pt_buffer<_locr,_locl>::Vector::iterator i = begin; i != end; ++i) {
		if(i != begin && (i->global_diagonal() - (i-1)->global_diagonal()) <= padding[frame]) {
			stat.inc(Statistics::DUPLICATES);
			continue;
		}
		local.push_back(local_match<_val> (i->seed_offset_, ref->data(i->subject_)));
		floating_sw(&query[i->seed_offset_],
				local.back(),
				padding[frame],
				score_matrix::get().rawscore(program_options::gapped_xdrop),
				program_options::gap_open + program_options::gap_extend,
				program_options::gap_extend,
				transcript_buf,
				Traceback ());
		const int score = local.back().score_;
		std::pair<size_t,size_t> l = ref_seqs<_val>::data_->local_position(i->subject_);
		matches.push_back(Segment<_val> (score, frame, &local.back(), l.first));
		anchored_transform(local.back(), l.second, i->seed_offset_);
		stat.inc(Statistics::ALIGNED_QLEN, local.back().query_len_);

		//local.back().print(query, ref_seqs<_val>::get()[l.first], transcript_buf);

		to_source_space(local.back(), frame, dna_len);
		stat.inc(Statistics::SCORE_TOTAL, local.back().score_);
		stat.inc(Statistics::OUT_HITS);
	}
}

#endif /* ALIGN_SEQUENCE_H_ */

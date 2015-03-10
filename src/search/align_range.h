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

#ifndef ALIGN_RANGE_H_
#define ALIGN_RANGE_H_

#include "align.h"
#include "../basic/statistics.h"

template<typename _val, typename _locr, typename _locq, typename _locl>
void align_range(_locq q_pos,
				 const typename sorted_list<_locr>::Type::const_iterator &s,
				 Statistics &stats,
				 typename Trace_pt_buffer<_locr,_locl>::Iterator &out,
				 unsigned sid)
{
	unsigned i = 0;

	const _val* query = query_seqs<_val>::data_->data(q_pos);
	hit_filter<_val,_locr,_locq,_locl> hf (stats, q_pos, out);

	if(s.n <= program_options::hit_cap) {
		stats.inc(Statistics::SEED_HITS, s.n);
		while(i < s.n) {
			align<_val,_locr,_locq,_locl>(q_pos, query, s[i], stats, sid, hf);
			++i;
		}
	} else {
		while(i < s.n && s[i] != 0) {
			assert(position_filter(s[i], filter_treshold(s.n), s.key()));
			align<_val,_locr,_locq,_locl>(q_pos, query, s[i], stats, sid, hf);
			stats.inc(Statistics::SEED_HITS);
			++i;
		}
	}

	hf.finish();
}

template<typename _val, typename _locr, typename _locq, typename _locl>
void align_range(const typename sorted_list<_locq>::Type::const_iterator &q,
				 const typename sorted_list<_locr>::Type::const_iterator &s,
				 Statistics &stats,
				 typename Trace_pt_buffer<_locr,_locl>::Iterator &out,
				 const unsigned sid)
{
	for(unsigned i=0;i<q.n; ++i)
		align_range<_val,_locr,_locq,_locl>(_locq(q[i]), s, stats, out, sid);
}

template<typename _val, typename _locr, typename _locq, typename _locl>
void align_partition(unsigned hp,
		Statistics &stats,
		unsigned sid,
		typename sorted_list<_locr>::Type::const_iterator i,
		typename sorted_list<_locq>::Type::const_iterator j)
{
	typename Trace_pt_buffer<_locr,_locl>::Iterator out (*Trace_pt_buffer<_locr,_locl>::instance);
	while(!i.at_end() && !j.at_end() && !exception_state()) {
		if(i.key() < j.key()) {
			++i;
		} else if(j.key() < i.key()) {
			++j;
		} else {
			align_range<_val,_locr,_locq,_locl>(j, i, stats, out, sid);
			++i;
			++j;
		}
	}
}

#endif /* ALIGN_RANGE_H_ */

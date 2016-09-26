/****
Copyright (c) 2014-2016, University of Tuebingen, Benjamin Buchfink
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/

#ifndef ALIGN_RANGE_H_
#define ALIGN_RANGE_H_

#include "filter_hit.h"
#include "../basic/statistics.h"

inline void align_range(Loc q_pos,
				 const sorted_list::const_iterator &s,
				 Statistics &stats,
				 Trace_pt_buffer::Iterator &out,
				 unsigned sid)
{
	unsigned i = 0, n=0;

	const Letter* query = query_seqs::data_->data(q_pos);
	hit_filter hf (stats, q_pos, out);

	if(s.n <= config.hit_cap) {
		stats.inc(Statistics::SEED_HITS, s.n);
		while(i < s.n) {
			align(q_pos, query, s[i], stats, sid, hf);
			++i;
		}
	} else {
		while(i < s.n && s[i] != 0) {
			assert(position_filter(s[i], filter_treshold((unsigned)s.n), s.key()));
			align(q_pos, query, s[i], stats, sid, hf);
			stats.inc(Statistics::SEED_HITS);
			++i;
			++n;
		}
	}
#ifdef EXTRA
	//if(n > 64)
		//printf("%u\n",n);
#endif

	hf.finish();
}

inline void align_range(const sorted_list::const_iterator &q,
				 const sorted_list::const_iterator &s,
				 Statistics &stats,
				 Trace_pt_buffer::Iterator &out,
				 const unsigned sid)
{
#ifdef EXTRA
	//if(q.n > 4096)
		//printf("%lu %lu\n",q.n,s.n);
#endif
	for(unsigned i=0;i<q.n; ++i)
		align_range(Loc(q[i]), s, stats, out, sid);
}

struct Stage1_hit
{
	Stage1_hit(unsigned q_ref, unsigned q_offset, unsigned s_ref, unsigned s_offset) :
		q(q_ref + q_offset),
		s(s_ref + s_offset)
	{}
	bool operator<(const Stage1_hit &rhs) const
	{
		return q < rhs.q;
	}
	struct Query
	{
		unsigned operator()(const Stage1_hit &x) const
		{
			return x.q;
		}
	};
	unsigned q, s;
};

void search_seed(const sorted_list::const_iterator &q,
	const sorted_list::const_iterator &s,
	Statistics &stats,
	Trace_pt_buffer::Iterator &out,
	const unsigned sid);

void stage2_search(const sorted_list::const_iterator &q,
	const sorted_list::const_iterator &s,
	const vector<Stage1_hit> &hits,
	Statistics &stats,
	Trace_pt_buffer::Iterator &out,
	const unsigned sid);

inline void align_partition(unsigned hp,
		Statistics &stats,
		unsigned sid,
		sorted_list::const_iterator i,
		sorted_list::const_iterator j,
		unsigned thread_id)
{
#ifndef SIMPLE_SEARCH
	if (hp > 0)
		return;
#endif
	Trace_pt_buffer::Iterator* out = new Trace_pt_buffer::Iterator (*Trace_pt_buffer::instance, thread_id);
	while(!i.at_end() && !j.at_end()) {
		if(i.key() < j.key()) {
			++i;
		} else if(j.key() < i.key()) {
			++j;
		} else {
			if (!config.slow_search) {
				//cout << "n=" << stats.data_[Statistics::SEED_HITS] << endl;
				/*if (stats.data_[Statistics::SEED_HITS] > 10000000000lu)
				break;*/
				if(config.old_freq || i[0] != 0)
					search_seed(j, i, stats, *out, sid);
			} else
				align_range(j, i, stats, *out, sid);
			++i;
			++j;
		}
	}
	delete out;
}

#endif /* ALIGN_RANGE_H_ */

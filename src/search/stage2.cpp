/****
Copyright (c) 2016, Benjamin Buchfink
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

#include "align_range.h"
#include "../util/map.h"
#include "../dp/dp.h"
#include "../data/queries.h"
#include "hit_filter.h"
#include "../data/reference.h"
#include "collision.h"

#ifdef __SSE2__

void search_query_offset(Loc q,
	const sorted_list::const_iterator &s,
	vector<Stage1_hit>::const_iterator hits,
	vector<Stage1_hit>::const_iterator hits_end,
	Statistics &stats,
	Trace_pt_buffer::Iterator &out,
	const unsigned sid)
{
	const Letter* query = query_seqs::data_->data(q);
	hit_filter hf(stats, q, out);

	for (vector<Stage1_hit>::const_iterator i = hits; i < hits_end; ++i) {
		const Loc s_pos = s[i->s];
		const Letter* subject = ref_seqs::data_->data(s_pos);
		/*cout << sequence(query, 32) << endl;
		cout << sequence(subject, 32) << endl;*/

		unsigned delta, len;
		int score;
		if ((score = stage2_ungapped(query, subject, sid, delta, len)) < config.min_ungapped_raw_score)
			continue;

		stats.inc(Statistics::TENTATIVE_MATCHES2);
		
		if (!is_primary_hit(query - delta, subject - delta, delta, sid, len))
			continue;

		stats.inc(Statistics::TENTATIVE_MATCHES3);
		hf.push(s_pos, score);
	}

	hf.finish();
}

#else

void search_query_offset(Loc q,
	const sorted_list::const_iterator &s,
	vector<Stage1_hit>::const_iterator hits,
	vector<Stage1_hit>::const_iterator hits_end,
	Statistics &stats,
	Trace_pt_buffer::Iterator &out,
	const unsigned sid)
{
	const Letter* query = query_seqs::data_->data(q);
	unsigned q_num_ = std::numeric_limits<unsigned>::max(), seed_offset_;

	for (vector<Stage1_hit>::const_iterator i = hits; i < hits_end; ++i) {
		const Loc s_pos = s[i->s];
		const Letter* subject = ref_seqs::data_->data(s_pos);

		unsigned delta, len;
		int score;
		if ((score = stage2_ungapped(query, subject, sid, delta, len)) < config.min_ungapped_raw_score)
			continue;

		stats.inc(Statistics::TENTATIVE_MATCHES2);

#ifndef NO_COLLISION_FILTER
		if (!is_primary_hit(query - delta, subject - delta, delta, sid, len))
			continue;
#endif

		stats.inc(Statistics::TENTATIVE_MATCHES3);

		if (score < config.min_hit_raw_score) {
			const sequence s = ref_seqs::data_->fixed_window_infix(s_pos + config.seed_anchor);
			unsigned left;
			sequence query(query_seqs::data_->window_infix(q + config.seed_anchor, left));
			score = smith_waterman(query, s, config.hit_band, left, config.gap_open + config.gap_extend, config.gap_extend);
		}
		if (score >= config.min_hit_raw_score) {
			if (q_num_ == std::numeric_limits<unsigned>::max()) {
				std::pair<size_t, size_t> l(query_seqs::data_->local_position(q));
				q_num_ = (unsigned)l.first;
				seed_offset_ = (unsigned)l.second;
			}
			assert(subject < ref_seqs::get().raw_len());
			out.push(hit(q_num_, s_pos, seed_offset_));
			stats.inc(Statistics::TENTATIVE_MATCHES4);
		}
	}
}

#endif

void stage2_search(const sorted_list::const_iterator &q,
	const sorted_list::const_iterator &s,
	const vector<Stage1_hit> &hits,
	Statistics &stats,
	Trace_pt_buffer::Iterator &out,
	const unsigned sid)
{
	typedef Map<vector<Stage1_hit>::const_iterator, Stage1_hit::Query> Map_t;
	Map_t map(hits.begin(), hits.end());
	for (Map_t::Iterator i = map.begin(); i.valid(); ++i)
		search_query_offset(q[i.begin()->q], s, i.begin(), i.end(), stats, out, sid);
}
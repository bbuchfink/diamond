/****
DIAMOND protein aligner
Copyright (C) 2013-2019 Benjamin Buchfink <buchfink@gmail.com>

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

#include "search.h"
#include "../util/map.h"
#include "../data/queries.h"
#include "hit_filter.h"
#include "../data/reference.h"
#include "collision.h"
#include "../dp/dp_matrix.h"
#include "../dp/ungapped.h"
#include "../util/sequence/sequence.h"
#include "../dp/ungapped_simd.h"

namespace Search {
namespace DISPATCH_ARCH {

#ifdef __SSE2__

template<typename _score> thread_local vector<score_vector<_score>> DP_matrix<_score>::scores_;
template<typename _score> thread_local vector<score_vector<_score>> DP_matrix<_score>::hgap_;

thread_local vector<sequence> hit_filter::subjects_;

void search_query_offset(Loc q,
	const Packed_loc *s,
	vector<Stage1_hit>::const_iterator hits,
	vector<Stage1_hit>::const_iterator hits_end,
	Statistics &stats,
	Trace_pt_buffer::Iterator &out,
	const unsigned sid)
{
	if (hits_end <= hits)
		return;
	const Letter* query = query_seqs::data_->data(q);
	hit_filter hf(stats, q, out);

	if (config.fast_stage2) {

		const Letter *subjects[16];
		int scores[16];

		const sequence query_clipped = Util::Sequence::clip(query - config.ungapped_window, config.ungapped_window*2, config.ungapped_window);
		const int window_left = int(query - query_clipped.data()), window = (int)query_clipped.length();
		
		for (vector<Stage1_hit>::const_iterator i = hits; i < hits_end; i += 16) {

			const size_t n = std::min(vector<Stage1_hit>::const_iterator::difference_type(16), hits_end - i);
			for (size_t j = 0; j < n; ++j)
				subjects[j] = ref_seqs::data_->data(s[(i + j)->s]) - window_left;
			DP::DISPATCH_ARCH::window_ungapped(query_clipped.data(), subjects, n, window, scores);
			
			for (size_t j = 0; j < n; ++j)
				if (scores[j] >= config.min_ungapped_raw_score) {
					stats.inc(Statistics::TENTATIVE_MATCHES2);
					const sequence subject_clipped = Util::Sequence::clip(subjects[j], window, window_left);
					const int delta = int(subject_clipped.data() - subjects[j]);
					assert(delta <= window_left);
					if (is_primary_hit(query_clipped.data() + delta, subject_clipped.data(), window_left - delta, sid, (unsigned)subject_clipped.length())) {
						stats.inc(Statistics::TENTATIVE_MATCHES3);
						hf.push(s[(i + j)->s], scores[j]);
					}
				}
		}

	}
	else {

		for (vector<Stage1_hit>::const_iterator i = hits; i < hits_end; ++i) {
			const Loc s_pos = s[i->s];
			const Letter* subject = ref_seqs::data_->data(s_pos);

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

	}

	hf.finish();
}

#else

void search_query_offset(Loc q,
	const Packed_loc *s,
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
			score = smith_waterman(query, s, config.hit_band, left, score_matrix.gap_open() + score_matrix.gap_extend(), score_matrix.gap_extend());
		}
		if (score >= config.min_hit_raw_score) {
			if (q_num_ == std::numeric_limits<unsigned>::max()) {
				std::pair<size_t, size_t> l(query_seqs::data_->local_position(q));
				q_num_ = (unsigned)l.first;
				seed_offset_ = (unsigned)l.second;
			}
			out.push(hit(q_num_, s_pos, seed_offset_));
			stats.inc(Statistics::TENTATIVE_MATCHES4);
		}
	}
}

#endif

void stage2(const Packed_loc *q,
	const Packed_loc *s,
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

}}
#include "align_range.h"
#include "../util/map.h"

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

		unsigned delta, len;
		int score;
		if ((score = xdrop_ungapped(query, subject, shapes.get_shape(sid).length_, delta, len)) < config.min_ungapped_raw_score)
			continue;

		stats.inc(Statistics::TENTATIVE_MATCHES2);
		
#ifndef NO_COLLISION_FILTER
		if (!is_primary_hit(query - delta, subject - delta, delta, sid, len))
			continue;
#endif

		stats.inc(Statistics::TENTATIVE_MATCHES3);
		hf.push(s_pos, score);
	}

	hf.finish();
}

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
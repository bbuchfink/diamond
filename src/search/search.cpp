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

#include "align_range.h"
#include "hit_filter.h"
#include "sse_dist.h"

Trace_pt_buffer* Trace_pt_buffer::instance;

TLS_PTR vector<Finger_print> *Seed_filter::vq_ptr, *Seed_filter::vs_ptr;
TLS_PTR vector<Stage1_hit> *Seed_filter::hits_ptr;

const unsigned tile_size[] = { 1024, 128 };

#define FAST_COMPARE2(q, s, stats, q_ref, s_ref, q_offset, s_offset, hits) if (q.match(s) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1)
#define FAST_COMPARE(q, s, stats, q_ref, s_ref, q_offset, s_offset, hits) if (q.match(s) >= config.min_identities) hits.push_back(Stage1_hit(q_ref, q_offset, s_ref, s_offset))

void query_register_search(vector<Finger_print>::const_iterator q,
	vector<Finger_print>::const_iterator s,
	vector<Finger_print>::const_iterator s_end,
	const Range_ref &ref,
	vector<Stage1_hit> &hits,
	Statistics &stats)
{
	const unsigned q_ref = unsigned(q - ref.q_begin);
	unsigned s_ref = unsigned(s - ref.s_begin);
	Finger_print q1 = *(q++), q2 = *(q++), q3 = *(q++), q4 = *(q++), q5 = *(q++), q6 = *q;
	const vector<Finger_print>::const_iterator end2 = s_end - (s_end - s) % 4;
	for (; s < end2; ) {
		Finger_print s1 = *(s++), s2 = *(s++), s3 = *(s++), s4 = *(s++);
		stats.inc(Statistics::SEED_HITS, 6 * 4);
		FAST_COMPARE(q1, s1, stats, q_ref, s_ref, 0, 0, hits);
		FAST_COMPARE(q2, s1, stats, q_ref, s_ref, 1, 0, hits);
		FAST_COMPARE(q3, s1, stats, q_ref, s_ref, 2, 0, hits);
		FAST_COMPARE(q4, s1, stats, q_ref, s_ref, 3, 0, hits);
		FAST_COMPARE(q5, s1, stats, q_ref, s_ref, 4, 0, hits);
		FAST_COMPARE(q6, s1, stats, q_ref, s_ref, 5, 0, hits);
		FAST_COMPARE(q1, s2, stats, q_ref, s_ref, 0, 1, hits);
		FAST_COMPARE(q2, s2, stats, q_ref, s_ref, 1, 1, hits);
		FAST_COMPARE(q3, s2, stats, q_ref, s_ref, 2, 1, hits);
		FAST_COMPARE(q4, s2, stats, q_ref, s_ref, 3, 1, hits);
		FAST_COMPARE(q5, s2, stats, q_ref, s_ref, 4, 1, hits);
		FAST_COMPARE(q6, s2, stats, q_ref, s_ref, 5, 1, hits);
		FAST_COMPARE(q1, s3, stats, q_ref, s_ref, 0, 2, hits);
		FAST_COMPARE(q2, s3, stats, q_ref, s_ref, 1, 2, hits);
		FAST_COMPARE(q3, s3, stats, q_ref, s_ref, 2, 2, hits);
		FAST_COMPARE(q4, s3, stats, q_ref, s_ref, 3, 2, hits);
		FAST_COMPARE(q5, s3, stats, q_ref, s_ref, 4, 2, hits);
		FAST_COMPARE(q6, s3, stats, q_ref, s_ref, 5, 2, hits);
		FAST_COMPARE(q1, s4, stats, q_ref, s_ref, 0, 3, hits);
		FAST_COMPARE(q2, s4, stats, q_ref, s_ref, 1, 3, hits);
		FAST_COMPARE(q3, s4, stats, q_ref, s_ref, 2, 3, hits);
		FAST_COMPARE(q4, s4, stats, q_ref, s_ref, 3, 3, hits);
		FAST_COMPARE(q5, s4, stats, q_ref, s_ref, 4, 3, hits);
		FAST_COMPARE(q6, s4, stats, q_ref, s_ref, 5, 3, hits);
		s_ref += 4;
	}
	for (; s < s_end; ++s) {
		stats.inc(Statistics::SEED_HITS, 6);
		FAST_COMPARE(q1, *s, stats, q_ref, s_ref, 0, 0, hits);
		FAST_COMPARE(q2, *s, stats, q_ref, s_ref, 1, 0, hits);
		FAST_COMPARE(q3, *s, stats, q_ref, s_ref, 2, 0, hits);
		FAST_COMPARE(q4, *s, stats, q_ref, s_ref, 3, 0, hits);
		FAST_COMPARE(q5, *s, stats, q_ref, s_ref, 4, 0, hits);
		FAST_COMPARE(q6, *s, stats, q_ref, s_ref, 5, 0, hits);
		++s_ref;
	}
}

void inner_search(vector<Finger_print>::const_iterator q,
	vector<Finger_print>::const_iterator q_end,
	vector<Finger_print>::const_iterator s,
	vector<Finger_print>::const_iterator s_end,
	const Range_ref &ref,
	vector<Stage1_hit> &hits,
	Statistics &stats)
{
	unsigned q_ref = unsigned(q - ref.q_begin);
	for (; q < q_end; ++q) {
		unsigned s_ref = unsigned(s - ref.s_begin);
		for (vector<Finger_print>::const_iterator s2 = s; s2 < s_end; ++s2) {
			stats.inc(Statistics::SEED_HITS);
			FAST_COMPARE((*q), *s2, stats, q_ref, s_ref, 0, 0, hits);
			++s_ref;
		}
		++q_ref;
	}
}

void Seed_filter::tiled_search(vector<Finger_print>::const_iterator q,
	vector<Finger_print>::const_iterator q_end,
	vector<Finger_print>::const_iterator s,
	vector<Finger_print>::const_iterator s_end,
	const Range_ref &ref,
	unsigned level)
{
	switch (level) {
	case 0:
	case 1:
		for (; q < q_end; q += std::min(q_end - q, (ptrdiff_t)tile_size[level]))
			for (vector<Finger_print>::const_iterator s2 = s; s2 < s_end; s2 += std::min(s_end - s2, (ptrdiff_t)tile_size[level]))
				tiled_search(q, q + std::min(q_end - q, (ptrdiff_t)tile_size[level]), s2, s2 + std::min(s_end - s2, (ptrdiff_t)tile_size[level]), ref, level+1);
		break;
	case 2:
		for (; q < q_end; q += std::min(q_end-q,(ptrdiff_t)6))
			if (q_end - q < 6)
				inner_search(q, q_end, s, s_end, ref, hits, stats);
			else
				query_register_search(q, s, s_end, ref, hits, stats);
	}	
}

void load_fps(const sorted_list::const_iterator &i, vector<Finger_print> &v, const Sequence_set &seqs)
{
	v.clear();
	v.reserve(i.n);
	for (unsigned j = 0; j < i.n; ++j)
		v.push_back(Finger_print(seqs.data(i[j])));
}

void Seed_filter::run(const sorted_list::const_iterator &q, const sorted_list::const_iterator &s)
{
	hits.clear();
	load_fps(q, vq, *query_seqs::data_);
	load_fps(s, vs, *ref_seqs::data_);
	tiled_search(vq.begin(), vq.end(), vs.begin(), vs.end(), Range_ref(vq.begin(), vs.begin()), 0);
	std::sort(hits.begin(), hits.end());
	stats.inc(Statistics::TENTATIVE_MATCHES1, hits.size());
	stage2_search(q, s, hits, stats, out, sid);
}
/****
DIAMOND protein aligner
Copyright (C) 2020 Max Planck Society for the Advancement of Science e.V.

Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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

#include <limits.h>
#include "search.h"
#include "../util/map.h"
#include "../data/queries.h"
#include "../data/reference.h"
#include "../dp/ungapped.h"
#include "../util/sequence/sequence.h"
#include "../dp/ungapped_simd.h"
#include "../dp/dp.h"
#include "left_most.h"
#include "../util/simd/vector.h"
#include "../dp/scan_diags.h"
#include "finger_print.h"
#include "../util/algo/all_vs_all.h"

namespace Search {
namespace DISPATCH_ARCH {

static bool gapped_filter(const LongScoreProfile& query_profile, const Letter *subject, int query_len) {
	const int l = (int)query_profile.length();
	int scores[64];
	DP::scan_diags64(query_profile, sequence(subject, l), -32, 0, l, scores);
	return score_matrix.evalue_norm(DP::diag_alignment(scores, 64), query_len) <= config.gapped_filter_evalue2;
}

static double ungapped_evalue(int query_len) {
	return query_len >= 60 ? config.ungapped_evalue : 1000.0;
}

void search_query_offset(Loc q,
	const Packed_loc* s,
	vector<Stage1_hit>::const_iterator hits,
	vector<Stage1_hit>::const_iterator hits_end,
	Statistics& stats,
	Trace_pt_buffer::Iterator& out,
	const unsigned sid,
	const Context &context)
{
	// thread_local LongScoreProfile query_profile;

	constexpr auto N = vector<Stage1_hit>::const_iterator::difference_type(::DISPATCH_ARCH::SIMD::Vector<int8_t>::CHANNELS);
	if (hits_end <= hits)
		return;
	const Letter* query = query_seqs::data_->data(q);

	const Letter* subjects[N];
	int scores[N];

	const sequence query_clipped = Util::Sequence::clip(query - config.ungapped_window, config.ungapped_window * 2, config.ungapped_window);
	const int window_left = int(query - query_clipped.data()), window = (int)query_clipped.length();
	unsigned query_id = UINT_MAX, seed_offset = UINT_MAX;
	std::pair<size_t, size_t> l = query_seqs::data_->local_position(q);
	query_id = (unsigned)l.first;
	seed_offset = (unsigned)l.second;
	const int query_len = query_seqs::data_->length(query_id);
	const double ungapped_cutoff = ungapped_evalue(query_len);
	// bool ps = false;
	
	for (vector<Stage1_hit>::const_iterator i = hits; i < hits_end; i += N) {

		const size_t n = std::min(N, hits_end - i);
		for (size_t j = 0; j < n; ++j)
			subjects[j] = ref_seqs::data_->data(s[(i + j)->s]) - window_left;
		DP::window_ungapped_best(query_clipped.data(), subjects, n, window, scores);

		for (size_t j = 0; j < n; ++j)
			if (score_matrix.evalue_norm(scores[j], query_len) <= ungapped_cutoff) {
				stats.inc(Statistics::TENTATIVE_MATCHES2);
				if(left_most_filter(query_clipped, subjects[j], window_left, shapes[sid].length_, context, sid == 0, sid)) {
					stats.inc(Statistics::TENTATIVE_MATCHES3);
					//if (config.gapped_filter_evalue2 == 0.0)
						out.push(hit(query_id, s[(i + j)->s], seed_offset));
					/*else {
						if (!ps) {
							query_profile.set(query_clipped);
							ps = true;
						}
						if (gapped_filter(query_profile, subjects[j], query_len)) {
							stats.inc(Statistics::TENTATIVE_MATCHES4);
							out.push(hit(query_id, s[(i + j)->s], seed_offset));
						}
					}*/
				}
			}
	}

}

void search_query_offset(Loc q,
	const Packed_loc* s,
	const uint32_t *hits,
	const uint32_t *hits_end,
	Statistics& stats,
	Trace_pt_buffer::Iterator& out,
	const unsigned sid,
	const Context& context)
{
	// thread_local LongScoreProfile query_profile;

	constexpr auto N = vector<Stage1_hit>::const_iterator::difference_type(::DISPATCH_ARCH::SIMD::Vector<int8_t>::CHANNELS);
	const Letter* query = query_seqs::data_->data(q);

	const Letter* subjects[N];
	int scores[N];

	const sequence query_clipped = Util::Sequence::clip(query - config.ungapped_window, config.ungapped_window * 2, config.ungapped_window);
	const int window_left = int(query - query_clipped.data()), window = (int)query_clipped.length();
	unsigned query_id = UINT_MAX, seed_offset = UINT_MAX;
	std::pair<size_t, size_t> l = query_seqs::data_->local_position(q);
	query_id = (unsigned)l.first;
	seed_offset = (unsigned)l.second;
	const int query_len = query_seqs::data_->length(query_id);
	const double ungapped_cutoff = ungapped_evalue(query_len);
	const int score_cutoff = context.cutoff_table(query_len);
	// bool ps = false;

	for (const uint32_t *i = hits; i < hits_end; i += N) {

		const size_t n = std::min(N, hits_end - i);
		for (size_t j = 0; j < n; ++j)
			subjects[j] = ref_seqs::data_->data(s[*(i + j)]) - window_left;
		DP::window_ungapped_best(query_clipped.data(), subjects, n, window, scores);

		for (size_t j = 0; j < n; ++j)
			if (score_matrix.evalue_norm(scores[j], query_len) <= ungapped_cutoff) {
			//if(scores[j] > score_cutoff) {
				stats.inc(Statistics::TENTATIVE_MATCHES2);
				if (left_most_filter(query_clipped, subjects[j], window_left, shapes[sid].length_, context, sid == 0, sid)) {
					stats.inc(Statistics::TENTATIVE_MATCHES3);
					//if (config.gapped_filter_evalue2 == 0.0)
					out.push(hit(query_id, s[*(i + j)], seed_offset));
					/*else {
						if (!ps) {
							query_profile.set(query_clipped);
							ps = true;
						}
						if (gapped_filter(query_profile, subjects[j], query_len)) {
							stats.inc(Statistics::TENTATIVE_MATCHES4);
							out.push(hit(query_id, s[(i + j)->s], seed_offset));
						}
					}*/
				}
			}
	}

}

void stage2(const Packed_loc *q,
	const Packed_loc *s,
	const vector<Stage1_hit> &hits,
	Statistics &stats,
	Trace_pt_buffer::Iterator &out,
	const unsigned sid,
	const Context &context)
{
	typedef Map<vector<Stage1_hit>::const_iterator, Stage1_hit::Query> Map_t;
	Map_t map(hits.begin(), hits.end());
	for (Map_t::Iterator i = map.begin(); i.valid(); ++i)
		search_query_offset(q[i.begin()->q], s, i.begin(), i.end(), stats, out, sid, context);
}

struct Stage2 {

	Stage2(const Packed_loc* q, const Packed_loc* s, Statistics& stat, Trace_pt_buffer::Iterator& out, unsigned sid, const Context& context):
		q(q),
		s(s),
		stat(stat),
		out(out),
		sid(sid),
		context(context)
	{}

	void operator()(const FlatArray<uint32_t> &hits, uint32_t query_begin, uint32_t subject_begin) {
		stat.inc(Statistics::TENTATIVE_MATCHES1, hits.data_size());
		const uint32_t query_count = (uint32_t)hits.size();
		const Packed_loc* q_begin = q + query_begin, * s_begin = s + subject_begin;
		for (uint32_t i = 0; i < query_count; ++i) {
			const uint32_t* r1 = hits.begin(i), * r2 = hits.end(i);
			if (r2 == r1)
				continue;
			search_query_offset(q_begin[i], s_begin, r1, r2, stat, out, sid, context);
		}
	}

	const Packed_loc* q, * s;
	Statistics& stat;
	Trace_pt_buffer::Iterator& out;
	const unsigned sid;
	const Context& context;

};

typedef vector<Finger_print, Util::Memory::AlignmentAllocator<Finger_print, 16>> Container;

static void load_fps(const Packed_loc* p, size_t n, Container& v, const Sequence_set& seqs)
{
	v.clear();
	v.reserve(n);
	const Packed_loc* end = p + n;
	for (; p < end; ++p)
		v.emplace_back(seqs.data(*p));
}

void stage1(const Packed_loc* q, size_t nq, const Packed_loc* s, size_t ns, Statistics& stats, Trace_pt_buffer::Iterator& out, const unsigned sid, const Context& context)
{
	thread_local Container vq, vs;
	stats.inc(Statistics::SEED_HITS, nq * ns);
	load_fps(q, nq, vq, *query_seqs::data_);
	load_fps(s, ns, vs, *ref_seqs::data_);
	Stage2 callback(q, s, stats, out, sid, context);
	all_vs_all(vq.data(), vq.size(), vs.data(), vs.size(), config.tile_size, callback);
}

}}
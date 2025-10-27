/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

#include <limits.h>
#include "search.h"
#include "dp/ungapped.h"
#include "util/sequence/sequence.h"
#include "dp/ungapped_simd.h"
#include "left_most.h"
#include "run/config.h"
#include "util/simd/vector.h"
#include "data/block/block.h"

using std::vector;
using std::pair;
using std::tie;
using std::numeric_limits;

namespace Search {
namespace DISPATCH_ARCH {

static const int SHORT_QUERY_LEN = 85;

static int ungapped_cutoff(int query_len, const WorkSet& context) {
#ifdef UNGAPPED_SPOUGE
	return query_len > config.short_query_max_len ? context.cutoff_table(query_len, 50) : context.short_query_ungapped_cutoff;
#else
	if (context.cfg.ungapped_evalue == 0.0)
		return 0;
	else if (query_len <= config.short_query_max_len)
		return context.context.short_query_ungapped_cutoff;
	else if (query_len <= SHORT_QUERY_LEN && align_mode.query_translated)
		return context.cfg.cutoff_table_short(query_len);
	else
		return context.cfg.cutoff_table(query_len);
#endif
}

static int ungapped_window(int query_len) {
	if (query_len <= SHORT_QUERY_LEN && align_mode.query_translated)
		return query_len;
	else
		return config.ungapped_window;
}

static pair<unsigned, Loc> query_data(const PackedLoc q, const SequenceSet& query_seqs) {
	pair<size_t, size_t> l = query_seqs.local_position((uint64_t)q); // lazy eval?
	return { (unsigned)l.first, (unsigned)l.second };
}
static pair<unsigned, Loc> query_data(const PackedLocId q, const SequenceSet& query_seqs) {
	return { q.block_id, Loc((int64_t)q.pos - query_seqs.position(q.block_id, 0)) };
}

template<typename SeedLoc>
static void search_query_offset(const SeedLoc& q,
	const SeedLoc* s,
	FlatArray<uint32_t>::DataConstIterator hits,
	FlatArray<uint32_t>::DataConstIterator hits_end,
	WorkSet& work_set)
{
	constexpr int N = ::DISPATCH_ARCH::SIMD::Vector<int8_t>::LANES;
	const SequenceSet& ref_seqs = work_set.cfg.target->seqs(), &query_seqs = work_set.cfg.query->seqs();
	const Letter* query = query_seqs.data(q);

	const Letter* subjects[N];
	int scores[N];
	std::fill(scores, scores + N, INT_MAX);

	unsigned query_id = UINT_MAX;
	Loc seed_offset = numeric_limits<Loc>::max();
	tie(query_id, seed_offset) = query_data(q, query_seqs);
	const int query_len = query_seqs.length(query_id);
	const int score_cutoff = ungapped_cutoff(query_len, work_set);
	const int window = ungapped_window(query_len);
	const Sequence query_clipped = Util::Seq::clip(query - window, window * 2, window);
	const int window_left = int(query - query_clipped.data()), window_clipped = (int)query_clipped.length();
	const unsigned sid = work_set.shape_id;
	const bool chunked = work_set.cfg.index_chunks > 1;
	const unsigned hamming_filter_id = work_set.cfg.hamming_filter_id;
	const bool self = config.self && work_set.cfg.current_ref_block == 0;
	int n = 0;
	size_t hit_count = 0;

	const int interval_mod = config.left_most_interval > 0 ? seed_offset % config.left_most_interval : window_left, interval_overhang = std::max(window_left - interval_mod, 0);

	for (FlatArray<uint32_t>::DataConstIterator i = hits; i < hits_end; i += n) {

		n = std::min(N, int(hits_end - i));
		for (int j = 0; j < n; ++j)
			subjects[j] = ref_seqs.data(s[*(i + j)]) - window_left;
		if (score_cutoff)
			DP::window_ungapped_best(query_clipped.data(), subjects, n, window_clipped, scores);

#ifdef WITH_GLOBAL_RANKING
		if (config.global_ranking_targets)
			for (int j = 0; j < n; ++j)
				if (scores[j] == UCHAR_MAX)
					scores[j] = ::ungapped_window(query_clipped.data(), subjects[j], window_clipped);
#endif

		for (int j = 0; j < n; ++j) {
			if (scores[j] > score_cutoff) {
				if (self && block_id(s[*(i + j)]) == query_id)
					continue;
#ifdef UNGAPPED_SPOUGE
				std::pair<size_t, size_t> l = ref_seqs::data_->local_position(s[*(i + j)]);
				if (scores[j] < context.cutoff_table(query_len, ref_seqs::data_->length(l.first)))
					continue;
#endif
				work_set.stats.inc(Statistics::TENTATIVE_MATCHES2);
				if (work_set.cfg.minimizer_window || work_set.cfg.sketch_size || config.lin_stage1 || work_set.cfg.lin_stage1_target
					|| left_most_filter(query_clipped + interval_overhang, subjects[j] + interval_overhang, window_left - interval_overhang, shapes[sid].length_, work_set.context, sid == 0, sid, score_cutoff, chunked, hamming_filter_id)) {
					work_set.stats.inc(Statistics::TENTATIVE_MATCHES3);
					if (config.global_ranking_targets)
						work_set.global_ranking_buffer->write({ query_id, (uint64_t)s[*(i + j)], seed_offset, (uint16_t)scores[j], block_id(s[*(i + j)]) });
					else {
						if (hit_count++ == 0)
							work_set.out->new_query(query_id, seed_offset);
#ifdef HIT_KEEP_TARGET_ID
						work_set.out->write(query_id, (uint64_t)s[*(i + j)], (uint16_t)scores[j], block_id(s[*(i + j)]));
#else
						work_set.out->write(query_id, (uint64_t)s[*(i + j)], (uint16_t)scores[j]);
#endif
					}
				}
			}
		}
	}
}

template<typename SeedLoc>
static void FLATTEN search_tile(
	const FlatArray<uint32_t> &hits,
	int32_t query_begin,
	int32_t subject_begin,
	const SeedLoc* q,
	const SeedLoc* s,
	WorkSet& work_set)
{
	work_set.stats.inc(Statistics::TENTATIVE_MATCHES1, hits.data_size());
	const uint32_t query_count = (uint32_t)hits.size();
	const SeedLoc* q_begin = q + query_begin, *s_begin = s + subject_begin;
	for (uint32_t i = 0; i < query_count; ++i) {
		FlatArray<uint32_t>::DataConstIterator r1 = hits.begin(i), r2 = hits.end(i);
		if (r2 == r1)
			continue;
		search_query_offset(q_begin[i], s_begin, r1, r2, work_set);
	}
}

}}
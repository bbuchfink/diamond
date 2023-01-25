/****
DIAMOND protein aligner
Copyright (C) 2020-2022 Max Planck Society for the Advancement of Science e.V.

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
#include "../util/text_buffer.h"
#include "../util/memory/alignment.h"
#include "../run/config.h"

using std::vector;
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

static void search_query_offset(const SeedLoc& q,
	const SeedLoc* s,
	FlatArray<uint32_t>::DataConstIterator hits,
	FlatArray<uint32_t>::DataConstIterator hits_end,
	WorkSet& work_set)
{
	constexpr int N = ::DISPATCH_ARCH::SIMD::Vector<int8_t>::CHANNELS;
	const SequenceSet& ref_seqs = work_set.cfg.target->seqs(), &query_seqs = work_set.cfg.query->seqs();
	const Letter* query = query_seqs.data(q);

	const Letter* subjects[N];
	int scores[N];
	std::fill(scores, scores + N, INT_MAX);

	unsigned query_id = UINT_MAX;
	Loc seed_offset = numeric_limits<Loc>::max();
#ifdef KEEP_TARGET_ID
	query_id = q.block_id;
	seed_offset = Loc((int64_t)q.pos - query_seqs.position(query_id, 0));
#else
	std::pair<size_t, size_t> l = query_seqs.local_position((uint64_t)q); // lazy eval?
	query_id = (unsigned)l.first;
	seed_offset = (unsigned)l.second;
#endif
	const int query_len = query_seqs.length(query_id);
	const int score_cutoff = ungapped_cutoff(query_len, work_set);
	const int window = ungapped_window(query_len);
	const Sequence query_clipped = Util::Seq::clip(query - window, window * 2, window);
	const int window_left = int(query - query_clipped.data()), window_clipped = (int)query_clipped.length();
	const unsigned sid = work_set.shape_id;
	const bool chunked = work_set.cfg.index_chunks > 1;
	const unsigned hamming_filter_id = work_set.cfg.hamming_filter_id;
	int n = 0;
	size_t hit_count = 0;

	const int interval_mod = config.left_most_interval > 0 ? seed_offset % config.left_most_interval : window_left, interval_overhang = std::max(window_left - interval_mod, 0);

	for (FlatArray<uint32_t>::DataConstIterator i = hits; i < hits_end; i += n) {

		n = std::min(N, int(hits_end - i));
		for (int j = 0; j < n; ++j)
			subjects[j] = ref_seqs.data(s[*(i + j)]) - window_left;
		if (score_cutoff)
			DP::window_ungapped_best(query_clipped.data(), subjects, n, window_clipped, scores);

		if (config.global_ranking_targets)
			for (int j = 0; j < n; ++j)
				if (scores[j] == UCHAR_MAX)
					scores[j] = ::ungapped_window(query_clipped.data(), subjects[j], window_clipped);

		for (int j = 0; j < n; ++j) {
			if (scores[j] > score_cutoff) {
#ifdef UNGAPPED_SPOUGE
				std::pair<size_t, size_t> l = ref_seqs::data_->local_position(s[*(i + j)]);
				if (scores[j] < context.cutoff_table(query_len, ref_seqs::data_->length(l.first)))
					continue;
#endif
				work_set.stats.inc(Statistics::TENTATIVE_MATCHES2);
				if (work_set.cfg.minimizer_window || config.sketch_size
					|| left_most_filter(query_clipped + interval_overhang, subjects[j] + interval_overhang, window_left - interval_overhang, shapes[sid].length_, work_set.context, sid == 0, sid, score_cutoff, chunked, hamming_filter_id)) {
					work_set.stats.inc(Statistics::TENTATIVE_MATCHES3);
					if (hit_count++ == 0)
						*work_set.out = SerializerTraits<Hit>::make_sentry(query_id, seed_offset);
#ifdef KEEP_TARGET_ID
					if(config.global_ranking_targets)
						//*work_set.out = { query_id, (uint64_t)s[*(i + j)].block_id, seed_offset, (uint16_t)scores[j] };
						*work_set.out = { query_id, (uint64_t)s[*(i + j)].pos, seed_offset, (uint16_t)scores[j], s[*(i + j)].block_id };
					else
#ifdef HIT_KEEP_TARGET_ID
						*work_set.out = { query_id, (uint64_t)s[*(i + j)].pos, seed_offset, (uint16_t)scores[j], s[*(i + j)].block_id };
#else
						*work_set.out = { query_id, (uint64_t)s[*(i + j)].pos, seed_offset, (uint16_t)scores[j] };
#endif
#else
					*work_set.out = { query_id, (uint64_t)s[*(i + j)], seed_offset, (uint16_t)scores[j] };
#endif
				}
			}
		}
	}
}

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

typedef vector<FingerPrint, Util::Memory::AlignmentAllocator<FingerPrint, 16>> Container;

static void all_vs_all(const FingerPrint* a, uint32_t na, const FingerPrint* b, uint32_t nb, FlatArray<uint32_t>& out, unsigned hamming_filter_id) {
	for (uint32_t i = 0; i < na; ++i) {
		const FingerPrint e = a[i];
		out.next();
		for (uint32_t j = 0; j < nb; ++j)
			if (e.match(b[j]) >= hamming_filter_id)
				out.push_back(j);
	}
}

static void all_vs_all_self(const FingerPrint* a, uint32_t na, FlatArray<uint32_t>& out, unsigned hamming_filter_id) {
	for (uint32_t i = 0; i < na; ++i) {
		const FingerPrint e = a[i];
		out.next();
		for (uint32_t j = i + 1; j < na; ++j)
			if (e.match(a[j]) >= hamming_filter_id)
				out.push_back(j);
	}
}

static void load_fps(const SeedLoc* p, size_t n, Container& v, const SequenceSet& seqs)
{
	v.clear();
	v.reserve(n);
	const SeedLoc* end = p + n;
	for (; p < end; ++p) {
		v.emplace_back(seqs.data(*p));
	}
}

void FLATTEN stage1(const SeedLoc* q, int32_t nq, const SeedLoc* s, int32_t ns, WorkSet& work_set)
{
#ifdef __APPLE__
	thread_local Container vq, vs;
#else
	Container& vq = work_set.vq, &vs = work_set.vs;
#endif
	
	const int32_t tile_size = config.tile_size;
	load_fps(s, ns, vs, work_set.cfg.target->seqs());
	work_set.stats.inc(Statistics::SEED_HITS, nq * ns);
	load_fps(q, nq, vq, work_set.cfg.query->seqs());
	const int32_t qs = (int32_t)vq.size(), ss = (int32_t)vs.size();
	for (int32_t i = 0; i < qs; i += tile_size) {
		for (int32_t j = 0; j < ss; j += tile_size) {
			work_set.hits.clear();
			all_vs_all(vq.data() + i, std::min(tile_size, qs - i), vs.data() + j, std::min(tile_size, ss - j), work_set.hits, work_set.cfg.hamming_filter_id);
			search_tile(work_set.hits, i, j, q, s, work_set);
		}
	}
}

void FLATTEN stage1_lin(const SeedLoc* q, int32_t nq, const SeedLoc* s, int32_t ns, WorkSet& work_set)
{
#ifdef __APPLE__
	thread_local Container vq, vs;
#else
	Container& vq = work_set.vq, &vs = work_set.vs;
#endif

	const int32_t tile_size = config.tile_size;
	load_fps(q, 1, vq, work_set.cfg.query->seqs());
	load_fps(s, ns, vs, work_set.cfg.target->seqs());
	work_set.stats.inc(Statistics::SEED_HITS, ns);
	
	const int32_t ss = (int32_t)vs.size();
	for (int32_t j = 0; j < ss; j += tile_size) {
		work_set.hits.clear();
		all_vs_all(vq.data(), 1, vs.data() + j, std::min(tile_size, ss - j), work_set.hits, work_set.cfg.hamming_filter_id);
		search_tile(work_set.hits, 0, j, q, s, work_set);
	}
}

void FLATTEN stage1_lin_ranked(const SeedLoc* q, int32_t nq, const SeedLoc* s, int32_t ns, WorkSet& work_set)
{
#ifdef __APPLE__
	thread_local Container vq, vs;
#else
	Container& vq = work_set.vq, &vs = work_set.vs;
#endif

	const int32_t tile_size = config.tile_size;
	const int32_t ranking = work_set.kmer_ranking->highest_ranking(q, q + nq);
	load_fps(q + ranking, 1, vq, work_set.cfg.query->seqs());
	load_fps(s, ns, vs, work_set.cfg.target->seqs());
	work_set.stats.inc(Statistics::SEED_HITS, ns);

	const int32_t ss = (int32_t)vs.size();
	for (int32_t j = 0; j < ss; j += tile_size) {
		work_set.hits.clear();
		all_vs_all(vq.data(), 1, vs.data() + j, std::min(tile_size, ss - j), work_set.hits, work_set.cfg.hamming_filter_id);
		search_tile(work_set.hits, ranking, j, q, s, work_set);
	}
}

void FLATTEN stage1_self(const SeedLoc* q, int32_t nq, const SeedLoc* s, int32_t ns, WorkSet& work_set)
{
#ifdef __APPLE__
	thread_local Container vs;
#else
	Container& vs = work_set.vs;
#endif

	const int32_t tile_size = config.tile_size;
	load_fps(s, ns, vs, work_set.cfg.target->seqs());

	work_set.stats.inc(Statistics::SEED_HITS, ns*(ns - 1) / 2);
	const int32_t ss = (int32_t)vs.size();
	for (int32_t i = 0; i < ss; i += tile_size) {
		work_set.hits.clear();
		all_vs_all_self(vs.data() + i, std::min(tile_size, ss - i), work_set.hits, work_set.cfg.hamming_filter_id);
		search_tile(work_set.hits, i, i, s, s, work_set);
		for (int32_t j = i + tile_size; j < ss; j += tile_size) {
			work_set.hits.clear();
			all_vs_all(vs.data() + i, std::min(tile_size, ss - i), vs.data() + j, std::min(tile_size, ss - j), work_set.hits, work_set.cfg.hamming_filter_id);
			search_tile(work_set.hits, i, j, s, s, work_set);
		}
	}
}

}}
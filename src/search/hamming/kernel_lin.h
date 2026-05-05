/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

#include "../stage2.h"
#include "data/block/block.h"
#include "kernel.h"
#include "data/sequence_file.h"

using std::vector;
using ::DISPATCH_ARCH::FingerPrint;

namespace Search { namespace DISPATCH_ARCH {

static void FLATTEN stage1_query_lin(const PackedLoc* __restrict q, uint_fast32_t nq, const PackedLoc* __restrict s, uint_fast32_t ns, WorkSet& work_set)
{
#ifdef __APPLE__
	thread_local Container vq, vs;
#else
	Container& vq = work_set.vq, &vs = work_set.vs;
#endif

	const uint_fast32_t tile_size = config.tile_size;
	::DISPATCH_ARCH::load_fps(q, 1, vq, work_set.cfg.query->seqs());
	::DISPATCH_ARCH::load_fps(s, ns, vs, work_set.cfg.target->seqs());
	work_set.stats.inc(Statistics::SEED_HITS, ns);
	
	const uint_fast32_t ss = (uint_fast32_t)vs.size();
	for (uint_fast32_t j = 0; j < ss; j += tile_size) {
		const uint_fast32_t ts = std::min(tile_size, ss - j);
		work_set.hits.init(1, ts);
		all_vs_all(vq.data(), 1, vs.data() + j, ts, work_set.hits, work_set.cfg.hamming_filter_id);
		search_tile(work_set.hits, 0, j, q, s, 0, work_set);
	}
}

static void FLATTEN stage1_query_lin_ranked(const PackedLocId* __restrict q, uint_fast32_t nq, const PackedLocId* __restrict s, uint_fast32_t ns, WorkSet& work_set)
{
#ifdef __APPLE__
	thread_local Container vq, vs;
#else
	Container& vq = work_set.vq, &vs = work_set.vs;
#endif

	const uint_fast32_t tile_size = config.tile_size;
	const bool rank_by_id = flag_any(work_set.cfg.db->flags(), SequenceFile::Flags::RANK_BY_SEQID);
	const uint_fast32_t ranking = work_set.kmer_ranking->highest_ranking(q, q + nq, rank_by_id ? work_set.cfg.query.get() : nullptr);
	::DISPATCH_ARCH::load_fps(q + ranking, 1, vq, work_set.cfg.query->seqs());
	::DISPATCH_ARCH::load_fps(s, ns, vs, work_set.cfg.target->seqs());
	work_set.stats.inc(Statistics::SEED_HITS, ns);

	const uint_fast32_t ss = (uint_fast32_t)vs.size();
	for (uint_fast32_t j = 0; j < ss; j += tile_size) {
		const uint_fast32_t ts = std::min(tile_size, ss - j);
		work_set.hits.init(1, ts);
		all_vs_all(vq.data(), 1, vs.data() + j, ts, work_set.hits, work_set.cfg.hamming_filter_id);
		search_tile(work_set.hits, ranking, j, q, s, 0, work_set);
	}
}

static void FLATTEN stage1_longest_combo_lin(const PackedLocId* __restrict q, uint_fast32_t nq, const PackedLocId* __restrict s, uint_fast32_t ns, WorkSet& work_set)
{
#ifdef __APPLE__
	thread_local Container vq, vs;
#else
	Container& vq = work_set.vq, & vs = work_set.vs;
#endif

	const uint_fast32_t tile_size = config.tile_size;
	const SequenceSet& queries = work_set.cfg.query->seqs(), &targets = work_set.cfg.target->seqs();
	bool query_pivot = true;
	uint_fast32_t pivot = 0;
	Loc pivot_len = queries.length(q[0].block_id);
	for (uint_fast32_t i = 1; i < nq; ++i) {
		const Loc len = queries.length(q[i].block_id);
		if (len > pivot_len) {
			pivot = i;
			pivot_len = len;
		}
	}
	for (uint_fast32_t i = 0; i < ns; ++i) {
		const Loc len = targets.length(s[i].block_id);
		if (len > pivot_len) {
			query_pivot = false;
			pivot = i;
			pivot_len = len;
		}
	}
	if (query_pivot) {
		::DISPATCH_ARCH::load_fps(q + pivot, 1, vq, queries);
		::DISPATCH_ARCH::load_fps(s, ns, vs, targets);
		work_set.stats.inc(Statistics::SEED_HITS, ns);

		const uint_fast32_t ss = (uint_fast32_t)vs.size();
		for (uint_fast32_t j = 0; j < ss; j += tile_size) {
			const uint_fast32_t ts = std::min(tile_size, ss - j);
			work_set.hits.init(1, ts);
			all_vs_all(vq.data(), 1, vs.data() + j, ts, work_set.hits, work_set.cfg.hamming_filter_id);
			search_tile(work_set.hits, pivot, j, q, s, 1, work_set);
		}
	}
	else {
		::DISPATCH_ARCH::load_fps(q, nq, vq, queries);
		::DISPATCH_ARCH::load_fps(s + pivot, 1, vs, targets);
		work_set.stats.inc(Statistics::SEED_HITS, nq);

		const uint_fast32_t qs = (uint_fast32_t)vq.size();
		for (uint_fast32_t j = 0; j < qs; j += tile_size) {
			const uint_fast32_t tq = std::min(tile_size, qs - j);
			work_set.hits.init(tq, 1);
			all_vs_all(vq.data() + j, tq, vs.data(), 1, work_set.hits, work_set.cfg.hamming_filter_id);
			search_tile(work_set.hits, j, pivot, q, s, 2, work_set);
		}
	}
}

template<typename SeedLoc>
static void FLATTEN stage1_target_lin(const SeedLoc* __restrict q, uint_fast32_t nq, const SeedLoc* __restrict s, uint_fast32_t ns, WorkSet& work_set)
{
#ifdef __APPLE__
	thread_local Container vq, vs;
#else
	Container& vq = work_set.vq, & vs = work_set.vs;
#endif

	const uint_fast32_t tile_size = config.tile_size;
	::DISPATCH_ARCH::load_fps(q, nq, vq, work_set.cfg.query->seqs());
	::DISPATCH_ARCH::load_fps(s, 1, vs, work_set.cfg.target->seqs());
	work_set.stats.inc(Statistics::SEED_HITS, nq);

	const uint_fast32_t qs = (uint_fast32_t)vq.size();
	for (uint_fast32_t j = 0; j < qs; j += tile_size) {
		const uint_fast32_t tq = std::min(tile_size, qs - j);
		work_set.hits.init(tq, 1);
		all_vs_all(vq.data() + j, tq, vs.data(), 1, work_set.hits, work_set.cfg.hamming_filter_id);
		search_tile(work_set.hits, j, 0, q, s, 0, work_set);
	}
}

static void FLATTEN stage1_mutual_cov_query_lin(const PackedLocId* __restrict q, uint_fast32_t nq, const PackedLocId* __restrict s, uint_fast32_t ns, WorkSet& work_set)
{
	const double mlr = work_set.cfg.min_length_ratio;
	const bool self = config.self && work_set.cfg.current_ref_block == 0;
#ifdef __APPLE__
	thread_local Container vq, vs;
#else
	Container& vq = work_set.vq, & vs = work_set.vs;
#endif

	::DISPATCH_ARCH::load_fps(s, ns, vs, work_set.cfg.target->seqs());
	::DISPATCH_ARCH::load_fps(q, nq, vq, work_set.cfg.query->seqs());

	const uint_fast32_t qs = (uint_fast32_t)vq.size(), ss = (uint_fast32_t)vs.size();
	uint_fast32_t j = 0;
	for (uint_fast32_t i = 0; i < qs;) {
		const Loc qlen = work_set.cfg.query->seqs().length(q[i].block_id);
		uint_fast32_t j1 = j;
		for (; j1 < ss; ++j1) {
			const Loc tlen = work_set.cfg.target->seqs().length(s[j1].block_id);
			if ((double)tlen / qlen < mlr)
				break;
		}
		const uint_fast32_t ts = j1 - j;
		if (ts > 0) {
			work_set.hits.init(1, ts);
			//all_vs_all(vq.data() + i, 1, vs.data() + j, j1 - j, work_set.hits, work_set.cfg.hamming_filter_id);
			const uint_fast32_t qpos = self ? i + (j1 - j) / 2 : i;
			all_vs_all(vq.data() + qpos, 1, vs.data() + j, ts, work_set.hits, work_set.cfg.hamming_filter_id);
			search_tile(work_set.hits, qpos, j, q, s, 0, work_set);
		}
		j = j1;
		if (j == ss)
			break;
		const Loc tlen = work_set.cfg.target->seqs().length(s[j].block_id);
		for (; i < qs; ++i) {
			const Loc qlen = work_set.cfg.query->seqs().length(q[i].block_id);
			if ((double)tlen / qlen >= mlr)
				break;
		}
	}
}

static void FLATTEN stage1_mutual_cov_target_lin(const PackedLocId* __restrict q, uint_fast32_t nq, const PackedLocId* __restrict s, uint_fast32_t ns, WorkSet& work_set)
{
	const double mlr = work_set.cfg.min_length_ratio;
#ifdef __APPLE__
	thread_local Container vq, vs;
#else
	Container& vq = work_set.vq, & vs = work_set.vs;
#endif

	::DISPATCH_ARCH::load_fps(s, ns, vs, work_set.cfg.target->seqs());
	::DISPATCH_ARCH::load_fps(q, nq, vq, work_set.cfg.query->seqs());

	const uint_fast32_t qs = (uint_fast32_t)vq.size(), ss = (uint_fast32_t)vs.size();
	uint_fast32_t i = 0;
	for (uint_fast32_t j = 0; j < ss;) {
		const Loc tlen = work_set.cfg.target->seqs().length(s[j].block_id);
		uint_fast32_t i1 = i;
		for (; i1 < qs; ++i1) {
			const Loc qlen = work_set.cfg.query->seqs().length(q[i1].block_id);
			if ((double)qlen / tlen < mlr)
				break;
		}
		const uint_fast32_t tq = i1 - i;
		work_set.hits.init(tq, 1);
		all_vs_all(vq.data() + i, tq, vs.data() + j, 1, work_set.hits, work_set.cfg.hamming_filter_id);
		search_tile(work_set.hits, i, j, q, s, 0, work_set);
		i = i1;
		if (i == qs)
			break;
		const Loc qlen = work_set.cfg.query->seqs().length(q[i].block_id);
		for (; j < ss; ++j) {
			const Loc tlen = work_set.cfg.target->seqs().length(s[j].block_id);
			if ((double)qlen / tlen >= mlr)
				break;
		}
	}
}

}}
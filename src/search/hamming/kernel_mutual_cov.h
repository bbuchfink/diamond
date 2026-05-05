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

#pragma once
#include "../stage2.h"
#include "data/block/block.h"

using std::array;

namespace Search { namespace DISPATCH_ARCH {
	
static void all_vs_all_mutual_cov(const PackedLocId* q, const PackedLocId* s, const array<char, 48>* __restrict a, uint32_t na, const array<char, 48>* __restrict b, uint32_t nb, HitField& out, unsigned hamming_filter_id, WorkSet& work_set) {
	uint32_t j0 = 0, j1 = 0;
	const double mlr = work_set.cfg.min_length_ratio;
	for (uint32_t i = 0; i < na; ++i) {
		const FingerPrint e(a[i]);
		const Loc qlen = work_set.cfg.query->seqs().length(q[i].block_id);
		while (j0 < nb) {
			const Loc tlen = work_set.cfg.target->seqs().length(s[j0].block_id);
			const double lr = (double)qlen / tlen;
			if (lr >= mlr)
				break;
			++j0;
		}
		j1 = std::max(j1, j0);
		for (uint32_t j = j0; j < j1; ++j) {
			out.set(i, j, e.match(FingerPrint(b[j])) >= hamming_filter_id);
		}
		for (; j1 < nb; ++j1) {
			const Loc tlen = work_set.cfg.target->seqs().length(s[j1].block_id);
			const double lr = (double)tlen / qlen;
			if (lr < mlr)
				break;
			out.set(i, j1, e.match(FingerPrint(b[j1])) >= hamming_filter_id);
		}
	}
}

static void all_vs_all_self_mutual_cov(const PackedLocId* q, const array<char, 48>* __restrict a, uint32_t na, HitField& out, unsigned hamming_filter_id, WorkSet& work_set) {
	const double mlr = work_set.cfg.min_length_ratio;
	for (uint32_t i = 0; i < na; ++i) {
		const FingerPrint e(a[i]);
		const Loc qlen = work_set.cfg.query->seqs().length(q[i].block_id);
		for (uint32_t j = i + 1; j < na; ++j) {
			const Loc tlen = work_set.cfg.query->seqs().length(q[j].block_id);
			if ((double)tlen / qlen < mlr)
				break;
			work_set.stats.inc(Statistics::SEED_HITS);
			out.set(i, j, e.match(FingerPrint(a[j])) >= hamming_filter_id);
		}
	}
}

static void FLATTEN stage1_mutual_cov(const PackedLocId* q, uint_fast32_t nq, const PackedLocId* s, uint_fast32_t ns, WorkSet& work_set)
{
#ifdef __APPLE__
	thread_local Container vq, vs;
#else
	Container& vq = work_set.vq, & vs = work_set.vs;
#endif

	const uint_fast32_t tile_size = config.tile_size;
	::DISPATCH_ARCH::load_fps(s, ns, vs, work_set.cfg.target->seqs());
	work_set.stats.inc(Statistics::SEED_HITS, nq * ns);
	::DISPATCH_ARCH::load_fps(q, nq, vq, work_set.cfg.query->seqs());
	const uint_fast32_t qs = (uint_fast32_t)vq.size(), ss = (uint_fast32_t)vs.size();
	for (uint_fast32_t i = 0; i < qs; i += tile_size) {
		for (uint_fast32_t j = 0; j < ss; j += tile_size) {
			const uint_fast32_t tq = std::min(tile_size, qs - i);
			const uint_fast32_t ts = std::min(tile_size, ss - j);
			work_set.hits.init(tq, ts);
			all_vs_all_mutual_cov(q + i, s + j, vq.data() + i, tq, vs.data() + j, ts, work_set.hits, work_set.cfg.hamming_filter_id, work_set);
			search_tile(work_set.hits, i, j, q, s, 0, work_set);
		}
	}
}

static void FLATTEN stage1_self_mutual_cov(const PackedLocId* q, uint_fast32_t nq, const PackedLocId* s, uint_fast32_t ns, WorkSet& work_set)
{
#ifdef __APPLE__
	thread_local Container vs;
#else
	Container& vs = work_set.vs;
#endif

	::DISPATCH_ARCH::load_fps(s, ns, vs, work_set.cfg.target->seqs());

	//work_set.stats.inc(Statistics::SEED_HITS, ns * (ns - 1) / 2);
	const uint_fast32_t ss = (uint_fast32_t)vs.size();
	work_set.hits.init(ss, ss);
	all_vs_all_self_mutual_cov(s, vs.data(), ss, work_set.hits, work_set.cfg.hamming_filter_id, work_set);
	search_tile(work_set.hits, 0, 0, s, s, 0, work_set);
}

}}
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

#include "../stage2.h"
#include "data/block/block.h"
#include "kernel.h"

using std::vector;

namespace Search { namespace DISPATCH_ARCH {

static void FLATTEN stage1_query_lin(const PackedLoc* __restrict q, int32_t nq, const PackedLoc* __restrict s, int32_t ns, WorkSet& work_set)
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

static void FLATTEN stage1_query_lin_ranked(const PackedLocId* __restrict q, int32_t nq, const PackedLocId* __restrict s, int32_t ns, WorkSet& work_set)
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

template<typename SeedLoc>
static void FLATTEN stage1_target_lin(const SeedLoc* __restrict q, int32_t nq, const SeedLoc* __restrict s, int32_t ns, WorkSet& work_set)
{
#ifdef __APPLE__
	thread_local Container vq, vs;
#else
	Container& vq = work_set.vq, & vs = work_set.vs;
#endif

	const int32_t tile_size = config.tile_size;
	load_fps(q, nq, vq, work_set.cfg.query->seqs());
	load_fps(s, 1, vs, work_set.cfg.target->seqs());
	work_set.stats.inc(Statistics::SEED_HITS, nq);

	const int32_t qs = (int32_t)vq.size();
	for (int32_t j = 0; j < qs; j += tile_size) {
		work_set.hits.clear();
		all_vs_all(vq.data() + j, std::min(tile_size, qs - j), vs.data(), 1, work_set.hits, work_set.cfg.hamming_filter_id);
		search_tile(work_set.hits, j, 0, q, s, work_set);
	}
}

static void FLATTEN stage1_mutual_cov_query_lin(const PackedLocId* __restrict q, int32_t nq, const PackedLocId* __restrict s, int32_t ns, WorkSet& work_set)
{
	const double mlr = work_set.cfg.min_length_ratio;
	const bool self = config.self && work_set.cfg.current_ref_block == 0;
#ifdef __APPLE__
	thread_local Container vq, vs;
#else
	Container& vq = work_set.vq, & vs = work_set.vs;
#endif

	load_fps(s, ns, vs, work_set.cfg.target->seqs());
	load_fps(q, nq, vq, work_set.cfg.query->seqs());

	const int32_t qs = (int32_t)vq.size(), ss = (int32_t)vs.size();
	int32_t j = 0;
	for (int32_t i = 0; i < qs;) {
		const Loc qlen = work_set.cfg.query->seqs().length(q[i].block_id);
		int32_t j1 = j;
		for (; j1 < ss; ++j1) {
			const Loc tlen = work_set.cfg.target->seqs().length(s[j1].block_id);
			if ((double)tlen / qlen < mlr)
				break;
		}
		work_set.hits.clear();
		//all_vs_all(vq.data() + i, 1, vs.data() + j, j1 - j, work_set.hits, work_set.cfg.hamming_filter_id);
		const int32_t qpos = self ? i + (j1 - j) / 2 : i;
		all_vs_all(vq.data() + qpos, 1, vs.data() + j, j1 - j, work_set.hits, work_set.cfg.hamming_filter_id);
		search_tile(work_set.hits, qpos, j, q, s, work_set);
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

static void FLATTEN stage1_mutual_cov_target_lin(const PackedLocId* __restrict q, int32_t nq, const PackedLocId* __restrict s, int32_t ns, WorkSet& work_set)
{
	const double mlr = work_set.cfg.min_length_ratio;
#ifdef __APPLE__
	thread_local Container vq, vs;
#else
	Container& vq = work_set.vq, & vs = work_set.vs;
#endif

	load_fps(s, ns, vs, work_set.cfg.target->seqs());
	load_fps(q, nq, vq, work_set.cfg.query->seqs());

	const int32_t qs = (int32_t)vq.size(), ss = (int32_t)vs.size();
	int32_t i = 0;
	for (int32_t j = 0; j < ss;) {
		const Loc tlen = work_set.cfg.target->seqs().length(s[j].block_id);
		int32_t i1 = i;
		for (; i1 < qs; ++i1) {
			const Loc qlen = work_set.cfg.query->seqs().length(q[i1].block_id);
			if ((double)qlen / tlen < mlr)
				break;
		}
		work_set.hits.clear();
		all_vs_all(vq.data() + i, i1 - i, vs.data() + j, 1, work_set.hits, work_set.cfg.hamming_filter_id);
		search_tile(work_set.hits, i, j, q, s, work_set);
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
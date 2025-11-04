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

#pragma once
#include "search.h"
#include "../stage2.h"
#include "data/block/block.h"
#include "kernel.h"

using std::vector;

namespace Search { namespace DISPATCH_ARCH {
	
static void all_vs_all_mutual_cov(const PackedLocId* q, const PackedLocId* s, const FingerPrint* a, uint32_t na, const FingerPrint* b, uint32_t nb, FlatArray<uint32_t>& out, unsigned hamming_filter_id, WorkSet& work_set) {
	uint32_t j0 = 0, j1 = 0;
	const double mlr = work_set.cfg.min_length_ratio;
	for (uint32_t i = 0; i < na; ++i) {
		const FingerPrint e = a[i];
		const Loc qlen = work_set.cfg.query->seqs().length(q[i].block_id);
		out.next();
		while (j0 < nb) {
			const Loc tlen = work_set.cfg.target->seqs().length(s[j0].block_id);
			const double lr = (double)qlen / tlen;
			if (lr >= mlr)
				break;
			++j0;
		}
		j1 = std::max(j1, j0);
		for (uint32_t j = j0; j < j1; ++j) {
			if (e.match(b[j]) >= hamming_filter_id)
				out.push_back(j);
		}
		for (; j1 < nb; ++j1) {
			const Loc tlen = work_set.cfg.target->seqs().length(s[j1].block_id);
			const double lr = (double)tlen / qlen;
			if (lr < mlr)
				break;
			if (e.match(b[j1]) >= hamming_filter_id)
				out.push_back(j1);
		}
	}
}

static void all_vs_all_self_mutual_cov(const PackedLocId* q, const FingerPrint* a, uint32_t na, FlatArray<uint32_t>& out, unsigned hamming_filter_id, WorkSet& work_set) {
	const double mlr = work_set.cfg.min_length_ratio;
	for (uint32_t i = 0; i < na; ++i) {
		const FingerPrint e = a[i];
		const Loc qlen = work_set.cfg.query->seqs().length(q[i].block_id);
		out.next();
		for (uint32_t j = i + 1; j < na; ++j) {
			const Loc tlen = work_set.cfg.query->seqs().length(q[j].block_id);
			if ((double)tlen / qlen < mlr)
				break;
			work_set.stats.inc(Statistics::SEED_HITS);
			if (e.match(a[j]) >= hamming_filter_id)
				out.push_back(j);
		}
	}
}

static void FLATTEN stage1_mutual_cov(const PackedLocId* q, int32_t nq, const PackedLocId* s, int32_t ns, WorkSet& work_set)
{
#ifdef __APPLE__
	thread_local Container vq, vs;
#else
	Container& vq = work_set.vq, & vs = work_set.vs;
#endif

	const int32_t tile_size = config.tile_size;
	load_fps(s, ns, vs, work_set.cfg.target->seqs());
	work_set.stats.inc(Statistics::SEED_HITS, nq * ns);
	load_fps(q, nq, vq, work_set.cfg.query->seqs());
	const int32_t qs = (int32_t)vq.size(), ss = (int32_t)vs.size();
	for (int32_t i = 0; i < qs; i += tile_size) {
		for (int32_t j = 0; j < ss; j += tile_size) {
			work_set.hits.clear();
			all_vs_all_mutual_cov(q + i, s + j, vq.data() + i, std::min(tile_size, qs - i), vs.data() + j, std::min(tile_size, ss - j), work_set.hits, work_set.cfg.hamming_filter_id, work_set);
			search_tile(work_set.hits, i, j, q, s, work_set);
		}
	}
}

static void FLATTEN stage1_self_mutual_cov(const PackedLocId* q, int32_t nq, const PackedLocId* s, int32_t ns, WorkSet& work_set)
{
#ifdef __APPLE__
	thread_local Container vs;
#else
	Container& vs = work_set.vs;
#endif

	load_fps(s, ns, vs, work_set.cfg.target->seqs());

	//work_set.stats.inc(Statistics::SEED_HITS, ns * (ns - 1) / 2);
	const int32_t ss = (int32_t)vs.size();
	work_set.hits.clear();
	all_vs_all_self_mutual_cov(s, vs.data(), ss, work_set.hits, work_set.cfg.hamming_filter_id, work_set);
	search_tile(work_set.hits, 0, 0, s, s, work_set);
}

}}
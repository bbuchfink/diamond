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
#include "../stage2.h"
#include "data/block/block.h"
#include "kernel.h"

using std::vector;

namespace Search { namespace DISPATCH_ARCH {
	
static void all_vs_all_self(const array<char, 48>* __restrict a, uint_fast32_t na, HitField& out, unsigned hamming_filter_id) {
	const uint_fast32_t na2 = na & ~uint_fast32_t(3);
	uint_fast32_t i = 0;
	for (; i < na2; i += 4) {
		const FingerPrint e1(a[i]);
		const FingerPrint e2(a[i + 1]);
		const FingerPrint e3(a[i + 2]);
		const FingerPrint e4(a[i + 3]);
		out.set(i, i + 1, e1.match(e2) >= hamming_filter_id);
		out.set(i, i + 2, e1.match(e3) >= hamming_filter_id);
		out.set(i, i + 3, e1.match(e4) >= hamming_filter_id);
		out.set(i + 1, i + 2, e2.match(e3) >= hamming_filter_id);
		out.set(i + 1, i + 3, e2.match(e4) >= hamming_filter_id);
		out.set(i + 2, i + 3, e3.match(e4) >= hamming_filter_id);
		for (uint_fast32_t j = i + 4; j < na; ++j) {
			const FingerPrint fa(a[j]);
			out.set(i, j, e1.match(fa) >= hamming_filter_id);
			out.set(i+1, j, e2.match(fa) >= hamming_filter_id);
			out.set(i+2, j, e3.match(fa) >= hamming_filter_id);
			out.set(i+3, j, e4.match(fa) >= hamming_filter_id);
		}
	}
	for (; i < na; ++i) {
		const FingerPrint e(a[i]);
		for (uint_fast32_t j = i + 1; j < na; ++j)
			out.set(i, j, e.match(FingerPrint(a[j])) >= hamming_filter_id);
	}
}

template<typename SeedLoc>
static void FLATTEN stage1_self(const SeedLoc* q, uint_fast32_t nq, const SeedLoc* s, uint_fast32_t ns, WorkSet& work_set)
{
#ifdef __APPLE__
	thread_local Container vs;
#else
	Container& vs = work_set.vs;
#endif

	const uint_fast32_t tile_size = config.tile_size;
	::DISPATCH_ARCH::load_fps(s, ns, vs, work_set.cfg.target->seqs());

	work_set.stats.inc(Statistics::SEED_HITS, ns * (ns - 1) / 2);
	const uint_fast32_t ss = (uint_fast32_t)vs.size();
	for (uint_fast32_t i = 0; i < ss; i += tile_size) {
		const uint_fast32_t tq = std::min(tile_size, ss - i);
		work_set.hits.init(tq, tq);
		all_vs_all_self(vs.data() + i, tq, work_set.hits, work_set.cfg.hamming_filter_id);
		search_tile(work_set.hits, i, i, s, s, work_set);
		for (uint_fast32_t j = i + tile_size; j < ss; j += tile_size) {
			const uint_fast32_t tq = std::min(tile_size, ss - i);
			const uint_fast32_t ts = std::min(tile_size, ss - j);
			work_set.hits.init(tq, ts);
			all_vs_all(vs.data() + i, tq, vs.data() + j, ts, work_set.hits, work_set.cfg.hamming_filter_id);
			search_tile(work_set.hits, i, j, s, s, work_set);
		}
	}
}

}}
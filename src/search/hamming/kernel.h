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
#include "search/search.h"
#include "finger_print.h"

using std::array;
using ::DISPATCH_ARCH::FingerPrint;

namespace Search { namespace DISPATCH_ARCH {

static void all_vs_all(const array<char, 48>* __restrict a, size_t na, const array<char, 48>* __restrict b, size_t nb, HitField& out, const unsigned hamming_filter_id) {
	const size_t na2 = na & ~size_t(3);
	size_t i = 0;
	for (; i < na2; i += 4) {
		const FingerPrint e1(a[i]);
		const FingerPrint e2(a[i + 1]);
		const FingerPrint e3(a[i + 2]);
		const FingerPrint e4(a[i + 3]);
		for (size_t j = 0; j < nb; ++j) {
			const FingerPrint fb(b[j]);
			out.set(i, j, e1.match(fb) >= hamming_filter_id);
			out.set(i+1, j, e2.match(fb) >= hamming_filter_id);
			out.set(i+2, j, e3.match(fb) >= hamming_filter_id);
			out.set(i+3, j, e4.match(fb) >= hamming_filter_id);
		}
	}
	for (; i < na; ++i) {
		const FingerPrint e(a[i]);
		for (size_t j = 0; j < nb; ++j)
			out.set(i, j, e.match(FingerPrint(b[j])) >= hamming_filter_id);
	}
}

template<typename SeedLoc>
static void FLATTEN stage1(const SeedLoc* __restrict q, int32_t nq, const SeedLoc* __restrict s, int32_t ns, WorkSet& work_set)
{
#ifdef __APPLE__
	thread_local Container vq, vs;
#else
	Container& vq = work_set.vq, & vs = work_set.vs;
#endif

	const int32_t tile_size = config.tile_size;
	::DISPATCH_ARCH::load_fps(s, ns, vs, work_set.cfg.target->seqs());
	work_set.stats.inc(Statistics::SEED_HITS, nq * ns);
	::DISPATCH_ARCH::load_fps(q, nq, vq, work_set.cfg.query->seqs());
	const int32_t qs = (int32_t)vq.size(), ss = (int32_t)vs.size();
	for (int32_t i = 0; i < qs; i += tile_size) {
		for (int32_t j = 0; j < ss; j += tile_size) {
			const size_t tq = std::min(tile_size, qs - i);
			const size_t ts = std::min(tile_size, ss - j);
			work_set.hits.init(tq, ts);
			all_vs_all(vq.data() + i, tq, vs.data() + j, ts, work_set.hits, work_set.cfg.hamming_filter_id);
			search_tile(work_set.hits, i, j, q, s, work_set);
		}
	}
}

}}
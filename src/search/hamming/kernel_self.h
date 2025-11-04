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
	
static void all_vs_all_self(const FingerPrint* a, uint32_t na, FlatArray<uint32_t>& out, unsigned hamming_filter_id) {
	for (uint32_t i = 0; i < na; ++i) {
		const FingerPrint e = a[i];
		out.next();
		for (uint32_t j = i + 1; j < na; ++j)
			if (e.match(a[j]) >= hamming_filter_id)
				out.push_back(j);
	}
}

template<typename SeedLoc>
static void FLATTEN stage1_self(const SeedLoc* q, int32_t nq, const SeedLoc* s, int32_t ns, WorkSet& work_set)
{
#ifdef __APPLE__
	thread_local Container vs;
#else
	Container& vs = work_set.vs;
#endif

	const int32_t tile_size = config.tile_size;
	load_fps(s, ns, vs, work_set.cfg.target->seqs());

	work_set.stats.inc(Statistics::SEED_HITS, ns * (ns - 1) / 2);
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

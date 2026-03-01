/****
Copyright (C) 2012-2026 Benjamin J. Buchfink <buchfink@gmail.com>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

	http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
****/
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include "search/search.h"
#include "finger_print.h"

using std::array;
using ::DISPATCH_ARCH::FingerPrint;

namespace Search { namespace DISPATCH_ARCH {

static void all_vs_all(const array<char, 48>* __restrict a, uint_fast32_t na, const array<char, 48>* __restrict b, uint_fast32_t nb, HitField& out, const unsigned hamming_filter_id) {
	const uint_fast32_t na2 = na & ~uint_fast32_t(3);
	uint_fast32_t i = 0;
	for (; i < na2; i += 4) {
		const FingerPrint e1(a[i]);
		const FingerPrint e2(a[i + 1]);
		const FingerPrint e3(a[i + 2]);
		const FingerPrint e4(a[i + 3]);
		for (uint_fast32_t j = 0; j < nb; ++j) {
			const FingerPrint fb(b[j]);
			out.set(i, j, e1.match(fb) >= hamming_filter_id);
			out.set(i+1, j, e2.match(fb) >= hamming_filter_id);
			out.set(i+2, j, e3.match(fb) >= hamming_filter_id);
			out.set(i+3, j, e4.match(fb) >= hamming_filter_id);
		}
	}
	for (; i < na; ++i) {
		const FingerPrint e(a[i]);
		for (uint_fast32_t j = 0; j < nb; ++j)
			out.set(i, j, e.match(FingerPrint(b[j])) >= hamming_filter_id);
	}
}

template<typename SeedLoc>
static void FLATTEN stage1(const SeedLoc* __restrict q, uint_fast32_t nq, const SeedLoc* __restrict s, uint_fast32_t ns, WorkSet& work_set)
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
			all_vs_all(vq.data() + i, tq, vs.data() + j, ts, work_set.hits, work_set.cfg.hamming_filter_id);
			search_tile(work_set.hits, i, j, q, s, work_set);
		}
	}
}

}}
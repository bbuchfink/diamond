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
#include <algorithm>
#include "search/search.h"
#include "finger_print.h"

using std::array;
using ::DISPATCH_ARCH::FingerPrint;

namespace Search { namespace DISPATCH_ARCH {

/* Compares every query fingerprint of a tile against every target fingerprint. The
   result bits of a query row are accumulated in a register 64 targets at a time and
   stored as whole words, so the inner loop carries no memory read-modify-write per pair. */
static void all_vs_all(const array<char, 48>* __restrict a, uint_fast32_t na, const array<char, 48>* __restrict b, uint_fast32_t nb, HitField& out, const unsigned hamming_filter_id) {
	const uint_fast32_t na2 = na & ~uint_fast32_t(3);
	uint_fast32_t i = 0;
	for (; i < na2; i += 4) {
		const FingerPrint e1(a[i]);
		const FingerPrint e2(a[i + 1]);
		const FingerPrint e3(a[i + 2]);
		const FingerPrint e4(a[i + 3]);
		uint64_t* const r1 = out.row(i), * const r2 = out.row(i + 1), * const r3 = out.row(i + 2), * const r4 = out.row(i + 3);
		for (uint_fast32_t j0 = 0; j0 < nb; j0 += 64) {
			const uint_fast32_t j1 = std::min(j0 + 64, nb);
			uint64_t m1 = 0, m2 = 0, m3 = 0, m4 = 0, bit = 1;
			for (uint_fast32_t j = j0; j < j1; ++j, bit <<= 1) {
				const FingerPrint fb(b[j]);
				m1 |= bit & -uint64_t(e1.match(fb) >= hamming_filter_id);
				m2 |= bit & -uint64_t(e2.match(fb) >= hamming_filter_id);
				m3 |= bit & -uint64_t(e3.match(fb) >= hamming_filter_id);
				m4 |= bit & -uint64_t(e4.match(fb) >= hamming_filter_id);
			}
			const size_t w = j0 >> 6;
			r1[w] = m1;
			r2[w] = m2;
			r3[w] = m3;
			r4[w] = m4;
		}
	}
	for (; i < na; ++i) {
		const FingerPrint e(a[i]);
		uint64_t* const r = out.row(i);
		for (uint_fast32_t j0 = 0; j0 < nb; j0 += 64) {
			const uint_fast32_t j1 = std::min(j0 + 64, nb);
			uint64_t m = 0, bit = 1;
			for (uint_fast32_t j = j0; j < j1; ++j, bit <<= 1)
				m |= bit & -uint64_t(e.match(FingerPrint(b[j])) >= hamming_filter_id);
			r[j0 >> 6] = m;
		}
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
			search_tile(work_set.hits, i, j, q, s, 0, work_set);
		}
	}
}

}}
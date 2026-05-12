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
#include <cstddef>
#include <algorithm>

template<typename It1, typename It2, typename Out, typename Cmp>
void batch_binary_search(It1 q0, It1 q1, It2 t0, It2 t1, Out out, Cmp cmp, ptrdiff_t ti = 0) {
	const ptrdiff_t nq = q1 - q0;
	if (nq == 0)
		return;
	if (nq == 1)
		*out++ = std::upper_bound(t0, t1, *q0, cmp) - t0 + ti - 1;

	const ptrdiff_t nt = t1 - t0;
	if (nt == 1) {
		for (It1 it = q0; it < q1; ++it)
			*out++ = ti;
		return;
	}
	const ptrdiff_t d = nt / 2;
	const It2 mid_t = t0 + d;
	const It1 mid_q = std::lower_bound(q0, q1, *mid_t, cmp);
	batch_binary_search(q0, mid_q, t0, mid_t, out, cmp, ti);
	batch_binary_search(mid_q, q1, mid_t, t1, out, cmp, ti + d);
}
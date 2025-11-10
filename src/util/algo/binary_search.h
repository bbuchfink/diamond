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
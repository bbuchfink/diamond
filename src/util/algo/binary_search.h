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
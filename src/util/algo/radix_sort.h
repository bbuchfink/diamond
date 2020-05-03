/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

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

#ifndef RADIX_SORT2_H_
#define RADIX_SORT2_H_

#include <algorithm>
#include "radix_cluster.h"
#include <vector>

/*template<typename _t>
void radix_sort(Relation<_t> &R, unsigned total_bits, _t *buf = 0)
{
	bool dealloc = false;
	if (!buf) {
		//buf = new _t[R.n];
		//buf = (_t*)new char[sizeof(_t)*R.n];
		buf = (_t*)MemoryPool::global().alloc(sizeof(_t)*R.n);
		dealloc = true;
	}
	unsigned *hst = new unsigned[1 << config.radix_bits];
	_t *in = R.data, *out = buf;
	for (int bits = (int)total_bits; bits > 0; bits -= config.radix_bits) {
		radix_cluster(Relation<_t>(in, R.n), total_bits - bits, out, hst);
		std::swap(in, out);
	}
	//if(dealloc) delete[] buf;
	if (dealloc) MemoryPool::global().free(buf);
	delete[] hst;
}*/




/*using std::vector;
using tthread::thread;

template<typename _t, typename _int_t, unsigned _radix = 8>
_int_t get_radix(const _t& x, unsigned shift)
{
	return ((_int_t)x >> shift) & (1 << _radix - 1);
}

template<typename _t, typename _int_t, unsigned _radix = 8>
void build_histogram(typename vector<_t>::const_iterator begin, typename vector<_t>::const_iterator end, unsigned* hst, unsigned shift)
{
	memset(hst, 0, sizeof(unsigned) * (1 << _radix));
	for (typename vector<_t>::const_iterator i = begin; i < end; ++i)
		++hst[get_radix(*i, shift)];
}

template<typename _t, typename _int_t, unsigned _radix = 8>
void scatter(typename vector<_t>::const_iterator begin, typename vector<_t>::const_iterator end, typename vector<typename vector<_t>::iterator>::iterator dst, unsigned shift)
{
	for (typename vector<_t>::const_iterator i = begin; i < end; ++i)
		*(dst[get_radix(*i, shift)]++) = *i;
}

template<typename _t, typename _int_t, unsigned _radix = 8>
void radix_sort(typename vector<_t>::iterator begin, typename vector<_t>::iterator end, unsigned n_threads)
{
	typedef unsigned Histogram[1 << _radix];
	const size_t n = end - begin;
	if (n == 0)
		return;
	vector<_t> buf(n);
	typename vector<_t>::iterator buf_begin = buf.begin();
	partition p(n, n_threads);

	typename vector<_t>::iterator* src = &begin, * dst = &buf_begin;

	for (unsigned shift = 0; shift < sizeof(_int_t) * 8 / _radix; shift += _radix) {
		vector<Histogram> hst(p.parts);
		vector<thread*> threads;
		for (unsigned i = 0; i < p.parts; ++i)
			threads.push_back(launch_thread(build_histogram, *src + p.getMin(i), *src + p.getMax(i), hst[i], shift));
		for (unsigned i = 0; i < p.parts; ++i) {
			threads[i]->join();
			delete threads[i];
		}
		vector<vector<typename vector<_t>::iterator> > pointers(p.parts);

		std::swap(src, dst);
	}
}*/

#endif
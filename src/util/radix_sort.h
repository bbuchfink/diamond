/****
Copyright (c) 2015, University of Tuebingen
Author: Benjamin Buchfink
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/

#ifndef RADIX_SORT_H_
#define RADIX_SORT_H_

#include <vector>
#include "thread.h"

using std::vector;
using tthread::thread;

template<typename _t, typename _int_t, unsigned _radix=8>
_int_t get_radix(const _t &x, unsigned shift)
{ return ((_int_t)x >> shift) & (1<<_radix-1); }

template<typename _t, typename _int_t, unsigned _radix=8>
void build_histogram(typename vector<_t>::const_iterator begin, typename vector<_t>::const_iterator end, unsigned *hst, unsigned shift)
{
	memset(hst, 0, sizeof(unsigned)*(1<<_radix));
	for(typename vector<_t>::const_iterator i=begin;i<end;++i)
		++hst[get_radix(*i, shift)];
}

template<typename _t, typename _int_t, unsigned _radix=8>
void scatter(typename vector<_t>::const_iterator begin, typename vector<_t>::const_iterator end, typename vector<typename vector<_t>::iterator>::iterator dst, unsigned shift)
{
	for(typename vector<_t>::const_iterator i=begin;i<end;++i)
		*(dst[get_radix(*i, shift)]++) = *i;
}

template<typename _t, typename _int_t, unsigned _radix=8>
void radix_sort(typename vector<_t>::iterator begin, typename vector<_t>::iterator end, unsigned n_threads)
{
	typedef unsigned Histogram[1<<_radix];
	const size_t n = end-begin;
	if(n == 0)
		return;
	vector<_t> buf (n);
	typename vector<_t>::iterator buf_begin = buf.begin();
	partition p (n, n_threads);

	typename vector<_t>::iterator *src = &begin, *dst = &buf_begin;

	for(unsigned shift=0; shift<sizeof(_int_t)*8/_radix; shift+=_radix) {
		vector<Histogram> hst (p.parts);
		vector<thread*> threads;
		for(unsigned i=0;i<p.parts;++i)
			threads.push_back(launch_thread(build_histogram, *src+p.getMin(i), *src+p.getMax(i), hst[i], shift));
		for(unsigned i=0;i<p.parts;++i) {
			threads[i]->join();
			delete threads[i];
		}
		vector<vector<typename vector<_t>::iterator> > pointers (p.parts);

		std::swap(src, dst);
	}
}

#endif /* RADIX_SORT_H_ */

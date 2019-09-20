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

template<typename _t>
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
}

#endif
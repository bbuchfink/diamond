/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#ifndef RADIX_CLUSTER_H_
#define RADIX_CLUSTER_H_

#include <string.h>

template<typename _t>
struct Relation
{
	Relation(_t *data, size_t n) :
		data(data),
		n(n)
	{}
	_t *data;
	size_t n;
	const _t& operator[](unsigned i) const
	{
		return data[i];
	}
	const _t *end() const
	{
		return data + n;
	}
};

struct ExtractBits
{
	ExtractBits(unsigned n, unsigned shift) :
		shift(shift),
		mask(n - 1)
	{}
	unsigned operator()(unsigned x) const
	{
		return (x >> shift) & mask;
	}
	unsigned shift, mask;
};

template<typename _t>
void radix_cluster(const Relation<_t> &in, unsigned shift, _t *&out, unsigned *&hst)
{
	static const size_t BUF_SIZE = 8;
	const unsigned clusters = 1 << config.radix_bits;
	ExtractBits radix(clusters, shift);

	hst = new unsigned[clusters];
	memset(hst, 0, clusters*sizeof(unsigned));
	for (const _t *i = in.data; i < in.end(); ++i)
		++hst[radix(i->key)];
	unsigned sum = 0;
	for (unsigned i = 0; i < clusters; ++i) {
		const unsigned c = hst[i];
		hst[i] = sum;
		sum += c;
	}

	out = new _t[in.n];

	if (config.radix_cluster_buffered) {
		_t *buf = new _t[clusters*BUF_SIZE];
		unsigned *buf_n = new unsigned[clusters];
		memset(buf_n, 0, clusters*sizeof(unsigned));

		for (const _t *i = in.data; i < in.end(); ++i) {
			const unsigned r = radix(i->key);
			buf[r*BUF_SIZE + buf_n[r]++] = *i;
			if (buf_n[r] == BUF_SIZE) {
				memcpy(out + hst[r], buf + r*BUF_SIZE, BUF_SIZE*sizeof(_t));
				hst[r] += BUF_SIZE;
				buf_n[r] = 0;
			}
		}
		delete[] buf;
		delete[] buf_n;
	}
	else {
		for (const _t *i = in.data; i < in.end(); ++i) {
			const unsigned r = radix(i->key);
			out[hst[r]++] = *i;
		}
	}
}

#endif
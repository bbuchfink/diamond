/****
DIAMOND protein aligner
Copyright (C) 2016-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
						
Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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

#pragma once
#include <thread>
#include "basic/config.h"
#include "partition.h"

template<typename T>
struct Relation
{
	Relation(T *data, size_t n) :
		data(data),
		n(n)
	{}
	T *data;
	size_t n;
	const T& operator[](size_t i) const
	{
		return data[i];
	}
	const T *end() const
	{
		return data + n;
	}
	Relation<T> part(size_t begin, size_t n) const {
		return Relation<T>(data + begin, n);
	}
};

template<typename T>
struct ExtractBits
{
	ExtractBits(T n, uint32_t shift) :
		shift(shift),
		mask(n - 1)
	{}
	T operator()(T x) const
	{
		return (x >> shift) & mask;
	}
	uint32_t shift;
	T mask;
};

template<typename T, typename GetKey>
void radix_cluster(const Relation<T> &in, unsigned shift, T *out, unsigned *hst)
{
	typedef typename T::Key Key;
	static const size_t BUF_SIZE = 8;
	const unsigned clusters = 1 << config.radix_bits;
	ExtractBits<Key> radix(clusters, shift);

	memset(hst, 0, clusters*sizeof(unsigned));
	for (const T *i = in.data; i < in.end(); ++i)
		++hst[radix(GetKey()(*i))];
	unsigned sum = 0;
	for (unsigned i = 0; i < clusters; ++i) {
		const unsigned c = hst[i];
		hst[i] = sum;
		sum += c;
	}

	if (config.radix_cluster_buffered) {
		T *buf = new T[clusters*BUF_SIZE];
		unsigned *buf_n = new unsigned[clusters];
		memset(buf_n, 0, clusters*sizeof(unsigned));

		for (const T *i = in.data; i < in.end(); ++i) {
			const unsigned r = radix(GetKey()(*i));
			buf[r*BUF_SIZE + buf_n[r]++] = *i;	
			if (buf_n[r] == BUF_SIZE) {
				memcpy(out + hst[r], buf + r*BUF_SIZE, BUF_SIZE*sizeof(T));
				hst[r] += BUF_SIZE;
				buf_n[r] = 0;
			}
		}
		delete[] buf;
		delete[] buf_n;
	}
	else {
		for (const T *i = in.data; i < in.end(); ++i) {
			const unsigned r = radix(GetKey()(*i));
			out[hst[r]++] = *i;
		}
	}
}

template<typename T, typename GetKey>
void parallel_radix_cluster_build_hst(Relation<T> in, uint32_t shift, size_t* hst)
{
	typedef typename T::Key Key;
	const Key clusters = (Key)1 << config.radix_bits;
	ExtractBits<Key> radix(clusters, shift);
	for (const T* i = in.data; i < in.end(); ++i)
		++hst[radix(GetKey()(*i))];
}

template<typename T, typename GetKey>
void parallel_radix_cluster_scatter(Relation<T> in, uint32_t shift, size_t* hst, T* out) {
	typedef typename T::Key Key;
	const Key clusters = (Key)1 << config.radix_bits;
	ExtractBits<Key> radix(clusters, shift);
	for (const T* i = in.data; i < in.end(); ++i) {
		const unsigned r = radix(GetKey()(*i));
		out[hst[r]++] = *i;
	}
}

template<typename T, typename GetKey>
void parallel_radix_cluster(const Relation<T>& in, uint32_t shift, T* out, size_t thread_count)
{
	typedef typename T::Key Key;
	const Key clusters = (Key)1 << config.radix_bits;
	std::vector<size_t> hst(clusters, 0);
		
	::Partition<size_t> p(in.n, thread_count);
	const size_t nt = p.parts;
	
	std::vector<std::vector<size_t>> thread_hst;
	thread_hst.reserve(nt);
	for (unsigned i = 0; i < nt; ++i)
		thread_hst.emplace_back(clusters, 0);
	std::vector<std::thread> threads;
	for (unsigned i = 0; i < nt; ++i) {
		threads.emplace_back(parallel_radix_cluster_build_hst<T, GetKey>, in.part(p.begin(i), p.size(i)), shift, thread_hst[i].data());
	}
	for (unsigned i = 0; i < nt; ++i) {
		threads[i].join();
		for (unsigned j = 0; j < clusters; ++j)
			hst[j] += thread_hst[i][j];
	}
	
	size_t sum = 0;
	for (unsigned i = 0; i < clusters; ++i) {
		const size_t c = hst[i];
		hst[i] = sum;
		sum += c;

		size_t sum2 = 0;
		for (size_t j = 0; j < nt; ++j) {
			const size_t d = thread_hst[j][i];
			thread_hst[j][i] = hst[i] + sum2;
			sum2 += d;
		}
	}

	threads.clear();
	for (unsigned i = 0; i < nt; ++i)
		threads.emplace_back(parallel_radix_cluster_scatter<T, GetKey>, in.part(p.begin(i), p.size(i)), shift, thread_hst[i].data(), out);
	for (unsigned i = 0; i < nt; ++i)
		threads[i].join();
}

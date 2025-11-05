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
#include <utility>
#include "basic/config.h"
#include "radix_cluster.h"
#include "../data_structures/hash_table.h"
#include "../data_structures/double_array.h"
#include "../math/integer.h"

struct RelPtr
{
	RelPtr()
	{}
	RelPtr(unsigned r):
		r(r),
		s(0)
	{}
	operator unsigned() const
	{
		return r;
	}
	unsigned r, s;
};

template<typename T>
void hash_table_join(
	const Relation<T> &R,
	const Relation<T> &S,
	unsigned shift,
	DoubleArray<typename T::Value> &dst_r,
	DoubleArray<typename T::Value> &dst_s)
{
	using Key = typename T::Key;
	typedef HashTable<Key, RelPtr, ExtractBits<Key>, NoModulo> Table;
	
	Key N = (Key)next_pow2(R.n * config.join_ht_factor);
	Table table(N, ExtractBits<Key>(N, shift));
	typename Table::Entry *p;

	for (T *i = R.data; i < R.end(); ++i) {
		p = table.insert(i->key);
		++p->value.r;
		i->key = Key(p - table.data());
	}

	T *hit_s = S.data;
	for (T *i = S.data; i < S.end(); ++i) {
		if ((p = table.find_entry(i->key))) {
			++p->value.s;
			hit_s->value = i->value;
			hit_s->key = Key(p - table.data());
			++hit_s;
		}
	}

	typename DoubleArray<typename T::Value>::Iterator it_r = dst_r.begin(), it_s = dst_s.begin();
	
	for (unsigned i = 0; i < table.size(); ++i) {
		p = &table.data()[i];
		if (p->value.s) {
			unsigned r = p->value.r, s = p->value.s;
			it_r.count() = r;
			it_s.count() = s;
			p->value.r = dst_r.offset(it_r) + 4;
			p->value.s = dst_s.offset(it_s) + 4;
			it_r.next();
			it_s.next();
		}
	}
	dst_r.set_end(it_r);
	dst_s.set_end(it_s);

	for (const T *i = R.data; i < R.end(); ++i) {
		p = &table.data()[i->key];
		if (p->value.s) {
			dst_r[p->value.r] = i->value;
			p->value.r += sizeof(typename T::Value);
		}
	}

	for (const T *i = S.data; i < hit_s; ++i) {
		p = &table.data()[i->key];
		dst_s[p->value.s] = i->value;
		p->value.s += sizeof(typename T::Value);
	}
}

template<typename T>
void table_join(
	const Relation<T> &R,
	const Relation<T> &S,
	unsigned total_bits,
	unsigned shift,
	DoubleArray<typename T::Value> &dst_r,
	DoubleArray<typename T::Value> &dst_s)
{
	using Key = typename T::Key;
	const Key keys = (Key)1 << (total_bits - shift);
	ExtractBits<Key> key(keys, shift);
	RelPtr *table = (RelPtr*)calloc(keys, sizeof(RelPtr));
	RelPtr *p;

	for (T *i = R.data; i < R.end(); ++i)
		++table[key(i->key)].r;

	T *hit_s = S.data;
	for (T *i = S.data; i < S.end(); ++i) {
		if ((p = &table[key(i->key)])->r) {
			++p->s;
			std::copy(i, i + 1, hit_s++);
		}
	}

	typename DoubleArray<typename T::Value>::Iterator it_r = dst_r.begin(), it_s = dst_s.begin();

	for (Key i = 0; i < keys; ++i) {
		p = &table[i];
		if (p->s) {
			unsigned r = p->r, s = p->s;
			it_r.count() = r;
			it_s.count() = s;
			p->r = dst_r.offset(it_r) + 4;
			p->s = dst_s.offset(it_s) + 4;
			it_r.next();
			it_s.next();
		}
	}
	dst_r.set_end(it_r);
	dst_s.set_end(it_s);

	for (const T *i = R.data; i < R.end(); ++i) {
		p = &table[key(i->key)];
		if (p->s) {
			dst_r[p->r] = i->value;
			p->r += sizeof(typename T::Value);
		}
	}

	for (const T *i = S.data; i < hit_s; ++i) {
		p = &table[key(i->key)];
		dst_s[p->s] = i->value;
		p->s += sizeof(typename T::Value);
	}

	free(table);
}

template<typename T>
void hash_join(
	Relation<T> R,
	Relation<T> S,
	T *dst_r,
	T *dst_s,
	DoubleArray<typename T::Value> &out_r,
	DoubleArray<typename T::Value> &out_s,
	unsigned total_bits = 32,
	unsigned shift = 0)
{
	if (R.n == 0 || S.n == 0)
		return;
	const unsigned key_bits = total_bits - shift;
	if (R.n < config.join_split_size || key_bits < config.join_split_key_len) {
		DoubleArray<typename T::Value> tmp_r((void*)dst_r), tmp_s((void*)dst_s);
		if (next_pow2(R.n * config.join_ht_factor) < 1llu << key_bits)
			hash_table_join(R, S, shift, tmp_r, tmp_s);
		else
			table_join(R, S, total_bits, shift, tmp_r, tmp_s);
		out_r.append(tmp_r);
		out_s.append(tmp_s);
	}
	else {
		const unsigned clusters = 1 << config.radix_bits;
		unsigned *hstR = new unsigned[clusters], *hstS = new unsigned[clusters];
		radix_cluster<T, typename T::GetKey>(R, shift, dst_r, hstR);
		radix_cluster<T, typename T::GetKey>(S, shift, dst_s, hstS);

		shift += config.radix_bits;
		hash_join(Relation<T>(dst_r, hstR[0]), Relation<T>(dst_s, hstS[0]), R.data, S.data, out_r, out_s, total_bits, shift);
		for (unsigned i = 1; i < clusters; ++i)
			hash_join(Relation<T>(dst_r + hstR[i - 1], hstR[i] - hstR[i - 1]), Relation<T>(dst_s + hstS[i - 1], hstS[i] - hstS[i - 1]), R.data + hstR[i - 1], S.data + hstS[i - 1], out_r, out_s, total_bits, shift);

		delete[] hstR;
		delete[] hstS;
	}
}

template<typename T>
std::pair<DoubleArray<typename T::Value>, DoubleArray<typename T::Value>> hash_join(Relation<T> R, Relation<T> S, unsigned total_bits = 32) {
	const bool swap = config.hash_join_swap && R.n > S.n;
	if (swap)
		std::swap(R, S);
	T *buf_r = (T*)malloc(sizeof(T) * R.n), *buf_s = (T*)malloc(sizeof(T) * S.n);
	DoubleArray<typename T::Value> out_r((void*)R.data), out_s((void*)S.data);
	hash_join(R, S, buf_r, buf_s, out_r, out_s, total_bits);
	free(buf_r);
	free(buf_s);
	if (swap)
		std::swap(out_r, out_s);
	return { out_r, out_s };
}

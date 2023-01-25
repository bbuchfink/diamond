/****
DIAMOND protein aligner
Copyright (C) 2016-2021 Max Planck Society for the Advancement of Science e.V.
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
#include <stdlib.h>
#include <cstdlib>
#include <utility>
#include <algorithm>
#include "../../basic/config.h"
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

template<typename _t>
void hash_table_join(
	const Relation<_t> &R,
	const Relation<_t> &S,
	unsigned shift,
	DoubleArray<typename _t::Value> &dst_r,
	DoubleArray<typename _t::Value> &dst_s)
{
	using Key = typename _t::Key;
	typedef HashTable<Key, RelPtr, ExtractBits<Key>, NoModulo> Table;
	
	Key N = (Key)next_power_of_2(R.n * config.join_ht_factor);
	Table table(N, ExtractBits<Key>(N, shift));
	typename Table::Entry *p;

	for (_t *i = R.data; i < R.end(); ++i) {
		p = table.insert(i->key);
		++p->value.r;
		i->key = Key(p - table.data());
	}

	_t *hit_s = S.data;
	for (_t *i = S.data; i < S.end(); ++i) {
		if ((p = table.find_entry(i->key))) {
			++p->value.s;
			hit_s->value = i->value;
			hit_s->key = Key(p - table.data());
			++hit_s;
		}
	}

	typename DoubleArray<typename _t::Value>::Iterator it_r = dst_r.begin(), it_s = dst_s.begin();
	
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

	for (const _t *i = R.data; i < R.end(); ++i) {
		p = &table.data()[i->key];
		if (p->value.s) {
			dst_r[p->value.r] = i->value;
			p->value.r += sizeof(typename _t::Value);
		}
	}

	for (const _t *i = S.data; i < hit_s; ++i) {
		p = &table.data()[i->key];
		dst_s[p->value.s] = i->value;
		p->value.s += sizeof(typename _t::Value);
	}
}

template<typename _t>
void table_join(
	const Relation<_t> &R,
	const Relation<_t> &S,
	unsigned total_bits,
	unsigned shift,
	DoubleArray<typename _t::Value> &dst_r,
	DoubleArray<typename _t::Value> &dst_s)
{
	using Key = typename _t::Key;
	const Key keys = (Key)1 << (total_bits - shift);
	ExtractBits<Key> key(keys, shift);
	RelPtr *table = (RelPtr*)calloc(keys, sizeof(RelPtr));
	RelPtr *p;

	for (_t *i = R.data; i < R.end(); ++i)
		++table[key(i->key)].r;

	_t *hit_s = S.data;
	for (_t *i = S.data; i < S.end(); ++i) {
		if ((p = &table[key(i->key)])->r) {
			++p->s;
			std::copy(i, i + 1, hit_s++);
		}
	}

	typename DoubleArray<typename _t::Value>::Iterator it_r = dst_r.begin(), it_s = dst_s.begin();

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

	for (const _t *i = R.data; i < R.end(); ++i) {
		p = &table[key(i->key)];
		if (p->s) {
			dst_r[p->r] = i->value;
			p->r += sizeof(typename _t::Value);
		}
	}

	for (const _t *i = S.data; i < hit_s; ++i) {
		p = &table[key(i->key)];
		dst_s[p->s] = i->value;
		p->s += sizeof(typename _t::Value);
	}

	free(table);
}

template<typename _t>
void hash_join(
	Relation<_t> R,
	Relation<_t> S,
	_t *dst_r,
	_t *dst_s,
	DoubleArray<typename _t::Value> &out_r,
	DoubleArray<typename _t::Value> &out_s,
	unsigned total_bits = 32,
	unsigned shift = 0)
{
	if (R.n == 0 || S.n == 0)
		return;
	const unsigned key_bits = total_bits - shift;
	if (R.n < config.join_split_size || key_bits < config.join_split_key_len) {
		DoubleArray<typename _t::Value> tmp_r((void*)dst_r), tmp_s((void*)dst_s);
		if (next_power_of_2(R.n * config.join_ht_factor) < 1llu << key_bits)
			hash_table_join(R, S, shift, tmp_r, tmp_s);
		else
			table_join(R, S, total_bits, shift, tmp_r, tmp_s);
		out_r.append(tmp_r);
		out_s.append(tmp_s);
	}
	else {
		const unsigned clusters = 1 << config.radix_bits;
		unsigned *hstR = new unsigned[clusters], *hstS = new unsigned[clusters];
		radix_cluster<_t, typename _t::GetKey>(R, shift, dst_r, hstR);
		radix_cluster<_t, typename _t::GetKey>(S, shift, dst_s, hstS);

		shift += config.radix_bits;
		hash_join(Relation<_t>(dst_r, hstR[0]), Relation<_t>(dst_s, hstS[0]), R.data, S.data, out_r, out_s, total_bits, shift);
		for (unsigned i = 1; i < clusters; ++i)
			hash_join(Relation<_t>(dst_r + hstR[i - 1], hstR[i] - hstR[i - 1]), Relation<_t>(dst_s + hstS[i - 1], hstS[i] - hstS[i - 1]), R.data + hstR[i - 1], S.data + hstS[i - 1], out_r, out_s, total_bits, shift);

		delete[] hstR;
		delete[] hstS;
	}
}

template<typename _t>
std::pair<DoubleArray<typename _t::Value>, DoubleArray<typename _t::Value>> hash_join(Relation<_t> R, Relation<_t> S, unsigned total_bits = 32) {
	const bool swap = config.hash_join_swap && R.n > S.n;
	if (swap)
		std::swap(R, S);
	_t *buf_r = (_t*)malloc(sizeof(_t) * R.n), *buf_s = (_t*)malloc(sizeof(_t) * S.n);
	DoubleArray<typename _t::Value> out_r((void*)R.data), out_s((void*)S.data);
	hash_join(R, S, buf_r, buf_s, out_r, out_s, total_bits);
	free(buf_r);
	free(buf_s);
	if (swap)
		std::swap(out_r, out_s);
	return { out_r, out_s };
}

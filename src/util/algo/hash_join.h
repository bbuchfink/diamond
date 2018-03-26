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

#ifndef HASH_JOIN_H_
#define HASH_JOIN_H_

#include <stdlib.h>
#include "../../basic/config.h"
#include "../util.h"
#include "radix_cluster.h"
#include "../data_structures/hash_table.h"
#include "../data_structures/double_array.h"

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
void hash_table_join(const Relation<_t> &R, const Relation<_t> &S, unsigned shift)
{
	typedef HashTable<unsigned, RelPtr, ExtractBits> Table;

	uint32_t N = next_power_of_2(R.n * config.join_ht_factor);

	Table table(N, ExtractBits(N, shift));
	typename Table::Entry *p;

	for (_t *i = R.data; i < R.end(); ++i) {
		p = table.insert(i->key());
		++p->r;
		i->key() = p - table.data();
	}

	unsigned keys_hit = 0;
	_t *hit_s = S.data;
	for (_t *i = S.data; i < S.end(); ++i) {
		if (p = table.find_entry(i->key())) {
			++p->s;
			hit_s->value() = i->value();
			hit_s->key() = p - table.data();
			++hit_s;
			if (p->s == 1)
				++keys_hit;
		}
	}

	unsigned sum_r = 0, sum_s = 1;
	DoubleArray<typename _t::Value> hits_r(keys_hit), hits_s(keys_hit);
	unsigned *limits_r = hits_r.limits(), *limits_s = hits_s.limits();
	
	for (unsigned i = 0; i < table.size(); ++i) {
		p = &table.data()[i];
		if (p->r && p->s) {
			unsigned r = p->r, s = p->s;
			p->r = sum_r;
			p->s = sum_s;
			*(limits_r++) = sum_r;
			*(limits_s++) = sum_s - 1;
			sum_r += r;
			sum_s += s;
		}
	}

	hits_r.init(sum_r);
	for (const _t *i = R.data; i < R.end(); ++i) {
		p = &table.data()[i->key()];
		if (p->s)
			hits_r.data()[p->r++] = i->value();
	}

	hits_s.init(sum_s - 1);
	for (const _t *i = S.data; i < hit_s; ++i) {
		p = &table.data()[i->key()];
		hits_s.data()[p->s++ - 1] = i->value();
	}
}


template<typename _t>
void table_join(const Relation<_t> &R, const Relation<_t> &S, unsigned total_bits, unsigned shift)
{
	const unsigned keys = 1 << (total_bits - shift);
	ExtractBits key(keys, shift);
	RelPtr *table = (RelPtr*)calloc(keys, sizeof(RelPtr));
	RelPtr *p;

	for (_t *i = R.data; i < R.end(); ++i)
		++table[key(i->key())].r;

	unsigned keys_hit = 0;
	_t *hit_s = S.data;
	for (_t *i = S.data; i < S.end(); ++i) {
		if ((p = &table[key(i->key())])->r) {
			++p->s;
			memcpy(hit_s++, i, sizeof(_t));
			if (p->s == 1)
				++keys_hit;
		}
	}

	unsigned sum_r = 0, sum_s = 1;
	DoubleArray<typename _t::Value> hits_r(keys_hit), hits_s(keys_hit);
	unsigned *limits_r = hits_r.limits(), *limits_s = hits_s.limits();

	for (unsigned i = 0; i < keys; ++i) {
		p = &table[i];
		if (p->r && p->s) {
			unsigned r = p->r, s = p->s;
			p->r = sum_r;
			p->s = sum_s;
			*(limits_r++) = sum_r;
			*(limits_s++) = sum_s - 1;
			sum_r += r;
			sum_s += s;
		}
	}

	hits_r.init(sum_r);
	for (const _t *i = R.data; i < R.end(); ++i) {
		p = &table[key(i->key())];
		if (p->s)
			hits_r.data()[p->r++] = i->value();
	}

	hits_s.init(sum_s - 1);
	for (const _t *i = S.data; i < hit_s; ++i) {
		p = &table[key(i->key())];
		hits_s.data()[p->s++ - 1] = i->value();
	}
	
	free(table);
}

template<typename _t>
void hash_join(const Relation<_t> &R, const Relation<_t> &S, unsigned total_bits = 32, unsigned shift = 0)
{
	if (R.n == 0 || S.n == 0)
		return;
	const unsigned key_bits = total_bits - shift;
	if (R.n < config.join_split_size || key_bits < config.join_split_key_len) {
		if (next_power_of_2(R.n * config.join_ht_factor) < 1 << key_bits)
			hash_table_join(R, S, shift);
		else
			table_join(R, S, total_bits, shift);			
	}
	else {
		_t *outR, *outS;
		unsigned *hstR, *hstS;
		radix_cluster(R, shift, outR, hstR);
		radix_cluster(S, shift, outS, hstS);

		shift += config.radix_bits;
		hash_join(Relation<_t>(outR, hstR[0]), Relation<_t>(outS, hstS[0]), total_bits, shift);
		for (unsigned i = 1; i < 1 << config.radix_bits; ++i)
			hash_join(Relation<_t>(&outR[hstR[i - 1]], hstR[i] - hstR[i - 1]), Relation<_t>(&outS[hstS[i - 1]], hstS[i] - hstS[i - 1]), total_bits, shift);

		delete[] hstR;
		delete[] hstS;
		delete[] outR;
		delete[] outS;
	}
}

#endif
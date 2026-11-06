/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <utility>
#include "basic/config.h"
#include "util/memory/mem_profile.h"
#include "radix_cluster.h"
#include "../data_structures/hash_table.h"
#include "../data_structures/double_array.h"
#include "../math/integer.h"

static inline uint32_t checked_double_array_count(size_t n) {
	if (n > std::numeric_limits<uint32_t>::max())
		throw std::overflow_error("Hash join bucket size exceeds DoubleArray count limit.");
	return static_cast<uint32_t>(n);
}

struct RelPtr
{
	RelPtr()
	{}
	RelPtr(size_t r):
		r(r),
		s(0)
	{}
	operator size_t() const
	{
		return r;
	}
	size_t r, s;
};

template<typename Entry>
struct HashJoinTableTlsBuf {
	Entry *ptr = nullptr;
	size_t cap = 0;
	~HashJoinTableTlsBuf() { std::free(ptr); }
	Entry* get_zeroed(size_t n) {
		if (n > cap) {
			std::free(ptr);
			cap = std::max(n, cap * 2);
			ptr = (Entry*)std::malloc(sizeof(Entry) * cap);
		}
		std::memset(ptr, 0, sizeof(Entry) * n);
		return ptr;
	}
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
	static thread_local HashJoinTableTlsBuf<typename Table::Entry> tls_ht;
	Table table((char*)tls_ht.get_zeroed(N), (size_t)N * sizeof(typename Table::Entry), ExtractBits<Key>(N, shift));
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
	
	for (size_t i = 0; i < table.size(); ++i) {
		p = &table.data()[i];
		if (p->value.s) {
			size_t r = p->value.r, s = p->value.s;
			it_r.count() = checked_double_array_count(r);
			it_s.count() = checked_double_array_count(s);
			p->value.r = dst_r.offset(it_r) + DoubleArray<typename T::Value>::header_size;
			p->value.s = dst_s.offset(it_s) + DoubleArray<typename T::Value>::header_size;
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

struct RelPtrTlsBuf {
	RelPtr *ptr = nullptr;
	size_t cap = 0;
	~RelPtrTlsBuf() { std::free(ptr); }
	RelPtr* get_zeroed(size_t n) {
		if (n > cap) {
			std::free(ptr);
			cap = std::max(n, cap * 2);
			ptr = (RelPtr*)std::malloc(sizeof(RelPtr) * cap);
		}
		std::memset(ptr, 0, sizeof(RelPtr) * n);
		return ptr;
	}
};

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
	static thread_local RelPtrTlsBuf tls_table;
	RelPtr *table = tls_table.get_zeroed(keys);
	RelPtr *p;

	const T *r_end = R.end();
	for (T *i = R.data; i < r_end; ++i)
		++table[key(i->key)].r;

	T *hit_s = S.data;
	const T *s_end = S.end();
	for (T *i = S.data; i < s_end; ++i) {
		if ((p = &table[key(i->key)])->r) {
			++p->s;
			std::copy(i, i + 1, hit_s++);
		}
	}

	typename DoubleArray<typename T::Value>::Iterator it_r = dst_r.begin(), it_s = dst_s.begin();

	for (Key i = 0; i < keys; ++i) {
		p = &table[i];
		if (p->s) {
			size_t r = p->r, s = p->s;
			it_r.count() = checked_double_array_count(r);
			it_s.count() = checked_double_array_count(s);
			p->r = dst_r.offset(it_r) + DoubleArray<typename T::Value>::header_size;
			p->s = dst_s.offset(it_s) + DoubleArray<typename T::Value>::header_size;
			it_r.next();
			it_s.next();
		}
	}
	dst_r.set_end(it_r);
	dst_s.set_end(it_s);

	for (const T *i = R.data; i < r_end; ++i) {
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
		size_t hstR_stack[256], hstS_stack[256];
		size_t *hstR = (clusters <= 256) ? hstR_stack : new size_t[clusters];
		size_t *hstS = (clusters <= 256) ? hstS_stack : new size_t[clusters];
		radix_cluster<T, typename T::GetKey>(R, shift, dst_r, hstR);
		radix_cluster<T, typename T::GetKey>(S, shift, dst_s, hstS);

		shift += config.radix_bits;
		hash_join(Relation<T>(dst_r, hstR[0]), Relation<T>(dst_s, hstS[0]), R.data, S.data, out_r, out_s, total_bits, shift);
		for (unsigned i = 1; i < clusters; ++i)
			hash_join(Relation<T>(dst_r + hstR[i - 1], hstR[i] - hstR[i - 1]), Relation<T>(dst_s + hstS[i - 1], hstS[i] - hstS[i - 1]), R.data + hstR[i - 1], S.data + hstS[i - 1], out_r, out_s, total_bits, shift);

		if (clusters > 256) {
			delete[] hstR;
			delete[] hstS;
		}
	}
}

template<typename T>
struct HashJoinScratchBuf {
	T *ptr = nullptr;
	size_t cap = 0;
	~HashJoinScratchBuf() { std::free(ptr); }
	T* get(size_t n) {
		if (n > cap) {
			std::free(ptr);
			cap = std::max(n, cap * 2);
			ptr = (T*)std::malloc(sizeof(T) * cap);
		}
		return ptr;
	}
};

template<typename T>
std::pair<DoubleArray<typename T::Value>, DoubleArray<typename T::Value>> hash_join(Relation<T> R, Relation<T> S, unsigned total_bits = 32) {
	const bool swap = config.hash_join_swap && R.n > S.n;
	if (swap)
		std::swap(R, S);
	static thread_local HashJoinScratchBuf<T> tls_buf_r, tls_buf_s;
	T *buf_r = tls_buf_r.get(R.n), *buf_s = tls_buf_s.get(S.n);
	DoubleArray<typename T::Value> out_r((void*)R.data), out_s((void*)S.data);
	hash_join(R, S, buf_r, buf_s, out_r, out_s, total_bits);
	if (swap)
		std::swap(out_r, out_s);
	return { out_r, out_s };
}

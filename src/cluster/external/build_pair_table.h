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
#include "util/util.h"
#include "util/algo/hash.h"
#include "external.h"
#include "input_buffer.h"

namespace External {

#pragma pack(1)

struct SeedEntry {
	static constexpr bool POD = true;
	SeedEntry() :
		seed(),
		oid(),
		len()
	{}
	SeedEntry(uint64_t seed, int64_t oid, int32_t len):
		seed(seed),
		oid(oid),
		len(len)
	{}
	uint64_t key() const {
		return hash64(seed);
	}
	bool operator<(const SeedEntry& e) const {
		return seed < e.seed || (seed == e.seed && (len > e.len || (len == e.len && oid < e.oid)));
	}
	friend void serialize(const SeedEntry& e, CompressedBuffer& buf) {
		buf.write((const char*)&e, sizeof(e));
	}
	friend void deserialize(InputFile& in, SeedEntry& e) {
		in.read(&e);
	}
	struct Key {
		uint64_t operator()(const SeedEntry& e) const {
			return e.seed;
		}
	};
	uint64_t seed;
	int64_t oid;
	int32_t len;
} PACKED_ATTRIBUTE;

#pragma pack()

inline void get_pairs_uni_cov(KeyMergeIterator<InputBuffer<SeedEntry>::ConstIterator, SeedEntry::Key>& it, BufferArray& buffers) {
	const int64_t rep_oid = it.begin()->oid;	
	const int32_t rep_len = it.begin()->len;
	for (auto j = it.begin() + 1; j < it.end(); ++j) {
		if (rep_oid == j->oid)
			continue;
		buffers.write_msb(PairEntry(rep_oid, j->oid, rep_len, j->len));
	}
}

inline void get_pairs_mutual_cov(KeyMergeIterator<InputBuffer<SeedEntry>::ConstIterator, SeedEntry::Key>& it, BufferArray& buffers) {
	const double mlr = config.min_length_ratio;
	const int32_t s = (int)it.count();
	int32_t j = 0;
	const auto begin = it.begin();
	for (int32_t i = 0; i < s;) {
		const Loc qlen = begin[i].len;
		int32_t j1 = j;
		for (; j1 < s; ++j1) {
			const Loc tlen = begin[j1].len;
			if ((double)tlen / qlen < mlr)
				break;
		}
		const int32_t qpos = i + (j1 - j) / 2;
		const int64_t rep_oid = begin[qpos].oid;
		const int32_t rep_len = begin[qpos].len;		
		for (int k = j; k < j1; ++k)
			if (rep_oid != begin[k].oid)
				buffers.write_msb(PairEntry(rep_oid, begin[k].oid, rep_len, begin[k].len));
		j = j1;
		if (j == s)
			break;
		const Loc tlen = begin[j].len;
		for (; i < s; ++i) {
			const Loc qlen = begin[i].len;
			if ((double)tlen / qlen >= mlr)
				break;
		}
	}
}

}
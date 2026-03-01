/****
Copyright (C) 2012-2026 Benjamin J. Buchfink <buchfink@gmail.com>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

	http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
****/
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include "util/util.h"
#include "util/algo/hash.h"
#include "external.h"
#include "input_buffer.h"

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
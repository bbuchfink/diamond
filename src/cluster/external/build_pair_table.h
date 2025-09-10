/****
Copyright © 2013-2025 Benjamin J. Buchfink

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
#include "util/util.h"
#include "util/hash_function.h"
#include "external.h"

namespace Cluster {

struct SeedEntry {
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
		return seed;
	}
	bool operator<(const SeedEntry& e) const {
		return seed < e.seed || (seed == e.seed && (len > e.len || (len == e.len && oid < e.oid)));
	}
	struct Key {
		uint64_t operator()(const SeedEntry& e) const {
			return e.seed;
		}
	};
	uint64_t seed;
	int64_t oid;
	int32_t len;
};

void get_pairs_uni_cov(KeyMergeIterator<const SeedEntry*, SeedEntry::Key>& it, BufferArray& buffers) {
	const int64_t rep_oid = it.begin()->oid, radix = MurmurHash()(rep_oid) & (RADIX_COUNT - 1);
	const int32_t rep_len = it.begin()->len;
	for (auto j = it.begin() + 1; j < it.end(); ++j) {
		if (rep_oid == j->oid)
			continue;
		buffers.write(radix, PairEntry(rep_oid, j->oid, rep_len, j->len));
	}
}

void get_pairs_mutual_cov(KeyMergeIterator<const SeedEntry*, SeedEntry::Key>& it, BufferArray& buffers) {
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
		const int64_t radix = MurmurHash()(rep_oid) & (RADIX_COUNT - 1);
		for (int k = j; k < j1; ++k)
			if (rep_oid != begin[k].oid)
				buffers.write(radix, PairEntry(rep_oid, begin[k].oid, rep_len, begin[k].len));
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
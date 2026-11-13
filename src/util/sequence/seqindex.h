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
#include <cstdint>
#include <vector>
#include "basic/sequence.h"
#include "basic/value.h"
#include "util/memory/memory_resource.h"
#include "basic/statistics.h"

struct Seqindex {

	Seqindex(Sequence seq, int k, double word_threshold, std::pmr::memory_resource& pool, Statistics& stats);	
	// Calls callback(query_pos, target_pos) for every k-mer of seq that matches a word of the index.
	template<typename F>
	void scan(Sequence seq, F&& callback) const {
		if (entry_.empty())
			return;
		const Letter* s = seq.data();
		const Loc len = seq.length();
		const Entry* const entries = entry_.data();
		const uint32_t* const head = head_.data();
		const uint16_t* const filter = filter_.data();
		uint32_t code = 0;
		int run = 0;
		for (Loc j = 0; j < len; ++j) {
			const Letter l = s[j];
			if (l < 0 || l >= TRUE_AA) {
				run = 0;
				code = 0;
				continue;
			}
			// Rolling code of the last k letters: the letter leaving the window is removed before shifting.
			if (run >= k_)
				code -= (uint32_t)s[j - k_] * lead_;
			code = code * TRUE_AA + (uint32_t)l;
			if (++run < k_)
				continue;
			const uint32_t h = hash(code), b = h >> shift_;
			// Most k-mers match no word, which the filter tells without a hardly predictable loop over the bucket.
			if ((filter[b] & filter_bit(h)) == 0)
				continue;
			const Loc begin = j - k_ + 1;
			for (const Entry* e = entries + head[b], *end = entries + head[b + 1]; e < end; ++e)
				if (e->code == code)
					callback(e->pos, begin);
		}
	}
	~Seqindex();
	int k() const {
		return k_;
	}

	static constexpr int MIN_KMER_LEN = 3;
	static constexpr int MAX_KMER_LEN = 6;

	struct Entry {
		uint32_t code;
		Loc pos;
	};

private:
	
	static uint32_t hash(uint32_t code) {
		return code * 2654435761u;
	}

	uint32_t bucket(uint32_t code) const {
		return hash(code) >> shift_;
	}

	// Bit of a hash value in the filter of its bucket, selected by the 4 bits following the bucket bits.
	uint16_t filter_bit(uint32_t h) const {
		return uint16_t(1u << ((h >> (shift_ - FILTER_BITS)) & ((1u << FILTER_BITS) - 1)));
	}

	static constexpr unsigned FILTER_BITS = 4;

	int k_;
	unsigned shift_;
	uint32_t lead_;                      // weight of the first letter of a k-mer in its code
	// The words grouped by hash bucket, in the order they were generated within a bucket. The words of bucket b
	// are entry_[head_[b]..head_[b + 1]), so that a lookup reads consecutive memory.
	std::pmr::vector<Entry> entry_;
	std::pmr::vector<uint32_t> head_;
	// Element b has the bits (see filter_bit) of the words of bucket b set.
	std::pmr::vector<uint16_t> filter_;

};
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
#include <vector>
#include "seed_histogram.h"
#include "search/seed_complexity.h"
#include "flags.h"

#ifndef __sparc__
#pragma pack(1)
#endif

struct Block;

template<typename SeedLoc>
struct SeedArray
{

	static PackedLoc make_seed_loc(PackedLoc pos, uint32_t block_id, PackedLoc) {
		return pos;
	}

	static PackedLocId make_seed_loc(PackedLoc pos, uint32_t block_id, PackedLocId) {
		return PackedLocId(pos, block_id);
	}

	struct Entry
	{
		Entry() :
			key(),
			value()
		{ }
		Entry(SeedOffset key, PackedLoc value) :
			key(key),
			value(value)
		{ }
		Entry(SeedOffset key, PackedLoc pos, uint32_t block_id) :
			key(key),
			value(make_seed_loc(pos, block_id, SeedLoc()))
		{}
		bool operator<(const Entry& entry)const {
			return  this->key < entry.key;
		}
		bool operator==(const Entry& entry)const {
			return this->key == entry.key && this->value == entry.value;
		}
		SeedOffset key;

		struct GetKey {
			uint32_t operator()(const Entry& e) const {
				return e.key;
			}
		};

		SeedLoc value;
		using Key = decltype(key);
		using Value = decltype(value);
		using value_type = Entry;
	} PACKED_ATTRIBUTE;


	template<typename Filter>
	SeedArray(Block &seqs, const ShapeHistogram &hst, const SeedPartitionRange &range, int seedp_bits, char *buffer, const Filter* filter, const EnumCfg& enum_cfg);

	template<typename Filter>
	SeedArray(Block& seqs, const SeedPartitionRange& range, int seedp_bits, const Filter* filter, EnumCfg& cfg);

	Entry* begin(unsigned i)
	{
		if (data_)
			return &data_[begin_[i]];
		else
			return entries_[i].data();
	}

	const Entry* begin(unsigned i) const
	{
		if (data_)
			return &data_[begin_[i]];
		else
			return entries_[i].data();
	}

	size_t size(size_t i) const
	{
		if (data_)
			return begin_[i + 1] - begin_[i];
		else {
			return entries_[i].size();
		}
	}

	size_t size() const {
		if (data_)
			return begin_.back();
		else {
			size_t n = 0;
			for (const auto& v : entries_)
				n += v.size();
			return n;
		}
	}

	const Search::SeedStats& stats() const {
		return stats_;
	}

	static char *alloc_buffer(const SeedHistogram &hst, int index_chunks);

	const int key_bits;

private:

	Entry *data_;
	std::vector<int64_t> begin_;
	std::vector<std::vector<Entry>> entries_;
	Search::SeedStats stats_;

};

#pragma pack()
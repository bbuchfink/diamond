/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

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
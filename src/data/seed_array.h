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
#include <array>
#include <vector>
#include "seed_histogram.h"
#include "../search/seed_complexity.h"
#include "flags.h"

#pragma pack(1)

struct Block;

struct SeedArray
{

	using Loc = PackedLoc;

	struct Entry
	{
		Entry() :
			key(),
			value()
		{ }
		Entry(SeedOffset key, Loc value) :
			key(key),
			value(value)
		{ }
		Entry(SeedOffset key, Loc pos, uint32_t block_id):
			key(key),
#ifdef KEEP_TARGET_ID
			value(pos, block_id)
#else
			value(pos)
#endif
		{}
		SeedOffset key;
		
		struct GetKey {
			uint32_t operator()(const Entry& e) const {
				return e.key;
			}
		};

		SeedLoc value;
		using Key = decltype(key);
		using Value = decltype(value);
	} PACKED_ATTRIBUTE;

	template<typename _filter>
	SeedArray(Block &seqs, const ShapeHistogram &hst, const SeedPartitionRange &range, char *buffer, const _filter *filter, const EnumCfg& enum_cfg);

	template<typename _filter>
	SeedArray(Block& seqs, const SeedPartitionRange& range, const _filter* filter, EnumCfg& cfg);

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
		if(data_)
			return begin_[Const::seedp];
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
	size_t begin_[Const::seedp + 1];
	std::array<std::vector<Entry>, Const::seedp> entries_;
	Search::SeedStats stats_;

};

#pragma pack()
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
#include "../basic/packed_loc.h"

#pragma pack(1)
// #define KEEP_TARGET_ID

struct SeedArray
{

	typedef PackedLoc _pos;

	struct Entry
	{

#ifdef KEEP_TARGET_ID
		struct Value {
			Value() {}
			Value(_pos pos) :
				pos(pos) {}
			Value(_pos pos, uint32_t block_id) :
				pos(pos),
				block_id(block_id)
			{}
			operator uint64_t() const {
				return (uint64_t)pos;
			}
			_pos pos;
			uint32_t block_id;
		};
#else
		typedef _pos Value;
#endif

		Entry() :
			key(),
			value()
		{ }
		Entry(unsigned key, _pos value) :
			key(key),
			value(value)
		{ }
		Entry(unsigned key, _pos pos, uint32_t block_id):
			key(key),
#ifdef KEEP_TARGET_ID
			value(pos, block_id)
#else
			value(pos)
#endif
		{}
		uint32_t key;
		
		struct GetKey {
			uint32_t operator()(const Entry& e) const {
				return e.key;
			}
		};

		typedef uint32_t Key;
		Value value;
	} PACKED_ATTRIBUTE;

	template<typename _filter>
	SeedArray(SequenceSet &seqs, size_t shape, const shape_histogram &hst, const SeedPartitionRange &range, const vector<size_t> &seq_partition, char *buffer, const _filter *filter, const SeedEncoding code, const std::vector<bool>* skip);

	template<typename _filter>
	SeedArray(SequenceSet& seqs, size_t shape, const SeedPartitionRange& range, const _filter* filter, const SeedEncoding code, const std::vector<bool>* skip);

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

	static char *alloc_buffer(const Partitioned_histogram &hst, size_t index_chunks);

	const size_t key_bits;

private:

	Entry *data_;
	size_t begin_[Const::seedp + 1];
	std::array<std::vector<Entry>, Const::seedp> entries_;

};

#pragma pack()
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

#ifndef SEED_ARRAY_H_
#define SEED_ARRAY_H_

#include "seed_histogram.h"
#include "../basic/packed_loc.h"

#pragma pack(1)

struct SeedArray
{

	typedef Packed_loc _pos;

	struct Entry
	{
		Entry() :
			key(),
			value()
		{ }
		Entry(unsigned key, _pos value) :
			key(key),
			value(value)
		{ }
		uint32_t key;
		_pos value;
		struct GetKey {
			uint32_t operator()(const Entry& e) const {
				return e.key;
			}
		};
		typedef _pos Value;
		typedef uint32_t Key;
	} PACKED_ATTRIBUTE;

	template<typename _filter>
	SeedArray(const Sequence_set &seqs, size_t shape, const shape_histogram &hst, const SeedPartitionRange &range, const vector<size_t> &seq_partition, char *buffer, const _filter *filter);

	Entry* begin(unsigned i)
	{
		return &data_[begin_[i]];
	}

	const Entry* begin(unsigned i) const
	{
		return &data_[begin_[i]];
	}

	size_t size(size_t i) const
	{
		return begin_[i + 1] - begin_[i];
	}

	static char *alloc_buffer(const Partitioned_histogram &hst);

private:

	Entry *data_;
	size_t begin_[Const::seedp + 1];

};

#pragma pack()

#endif
/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
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
		unsigned key;
		_pos value;
		typedef _pos Value;
	} PACKED_ATTRIBUTE;

	template<typename _filter>
	SeedArray(const Sequence_set &seqs, size_t shape, const shape_histogram &hst, const SeedPartitionRange &range, const vector<size_t> seq_partition, const _filter *filter);

	Entry* begin(unsigned i)
	{
		return data_[i];
	}

	const Entry* begin(unsigned i) const
	{
		return data_[i];
	}

	size_t size(size_t i) const
	{
		return size_[i];
	}

private:

	Entry *data_[Const::seedp];
	size_t size_[Const::seedp];

};

#pragma pack()

#endif
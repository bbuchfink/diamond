/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

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
#include "sequence_set.h"
#include "../util/hash_table.h"
#include "../util/ptr_vector.h"

struct Seed_set
{
	Seed_set(const Sequence_set &seqs, double max_coverage);
	bool contains(uint64_t key, uint64_t shape) const
	{
		return data_[key];
	}
	double coverage() const
	{
		return coverage_;
	}
private:
	std::vector<bool> data_;
	double coverage_;
};

struct Hashed_seed_set
{
	Hashed_seed_set(const Sequence_set &seqs);
	bool contains(uint64_t key, uint64_t shape) const
	{
		return data_[shape].contains(key);
	}
private:
	PtrVector<PHash_set<void, No_hash>> data_;
};

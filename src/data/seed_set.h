/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

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

#ifndef SEED_SET_H_
#define SEED_SET_H_

#include "sequence_set.h"
#include "../util/hash_table.h"

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
	vector<bool> data_;
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
	Ptr_vector<PHash_set<Modulo2, No_hash> > data_;
};

#endif
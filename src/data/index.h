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
#include <memory>
#include "sequence_set.h"
#include "../util/hash_table.h"
#include "seed_histogram.h"
#include "sorted_list.h"

using std::vector;

struct Seed_index
{

	Seed_index()
	{}
	Seed_index(const Partitioned_histogram &hst, const Sequence_set &seqs, unsigned sh);

	sorted_list::Random_access_iterator operator[](uint64_t seed) const
	{
		const unsigned p = seed_partition(seed);
		PHash_table<Entry>::entry *e = tables[p][murmur_hash()(seed_partition_offset(seed))];
		if (e == 0)
			return sorted_list::Random_access_iterator(0, 0);
		else
			return list.random_access(p, e->value.n);
	}
	
private:

	static const double hash_table_factor;
	
	struct Entry
	{
		uint32_t n;
		operator unsigned() const
		{
			return n;
		}
	};

	static void fill_tables(Atomic<unsigned> *seedp, Seed_index *idx);
	static void count_seeds(Atomic<unsigned> *seedp, vector<size_t> *counts, sorted_list *list);

	std::unique_ptr<char> list_buffer;
	sorted_list list;
	PHash_table<Entry> tables[Const::seedp];

};

extern Seed_index seed_index[Const::max_shapes];

vector<Array<unsigned, Hashed_seed::p> > count_exact(const Sequence_set &seqs);
vector<size_t> count_approximate(const Sequence_set &seqs);
void build_index(const Sequence_set &seqs);


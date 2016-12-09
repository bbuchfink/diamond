/****
Copyright (c) 2016, Benjamin Buchfink
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/

#ifndef INDEX_H_
#define INDEX_H_

#include <vector>
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

	auto_ptr<char> list_buffer;
	sorted_list list;
	PHash_table<Entry> tables[Const::seedp];

};

extern Seed_index seed_index[Const::max_shapes];

vector<Array<unsigned, Hashed_seed::p> > count_exact(const Sequence_set &seqs);
vector<size_t> count_approximate(const Sequence_set &seqs);
void build_index(const Sequence_set &seqs);

#endif
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

#include <numeric>
#include "index.h"
#include "queries.h"
#include "../util/log_stream.h"

const double Seed_double_index::hash_table_factor = 2;
Seed_double_index seed_index[Const::max_shapes];

Seed_double_index::Seed_double_index(size_t psize)
{
	for (unsigned p = 0; p < Hashed_seed::p; ++p)
		assign_ptr(tables[p], new PHash_table<Entry>(std::max(size_t((double)psize * hash_table_factor), (size_t)1llu)));
}

Seed_double_index::Seed_double_index(const Array<unsigned, Hashed_seed::p> &psize)
{
	for (unsigned p = 0; p < Hashed_seed::p; ++p)
		assign_ptr(tables[p], new PHash_table<Entry>(std::max(size_t((double)psize[p] * hash_table_factor), (size_t)1llu)));
}

struct Seed_entry
{
	Seed_entry()
	{}
	Seed_entry(Hashed_seed seed, size_t pos) :
		seed(seed)
	{}
	Hashed_seed seed;
};

struct Count_query_callback
{
	void operator()(unsigned shape_id, const Seed_entry& seed) const
	{
		PHash_table<Seed_double_index::Entry>::entry *e = seed_index[shape_id].tables[seed.seed.partition()].insert(seed.seed.offset());
		++e->value.q;
	}
};

struct Access
{
	void operator()(Hashed_seed seed, size_t pos, unsigned shape_id)
	{
		PHash_table<Seed_double_index::Entry>::entry *e = seed_index[shape_id].tables[seed.partition()][seed.offset()];
		if (e)
			++n;
	}
	void finish()
	{}
	unsigned n;
};

void build_query_index()
{
	static const size_t count_exact_limit = 10000000llu;

	const Sequence_set& seqs = query_seqs::get();

	task_timer timer("Counting query seeds", 2);
	vector<Array<unsigned, Hashed_seed::p> > exact_counts;
	vector<size_t> apxt_counts;
	const bool exact = seqs.letters() < count_exact_limit;
	if (exact) {
		exact_counts = count_exact(seqs);
		timer.finish();
		verbose_stream << "Seeds = " << std::accumulate(exact_counts[0].begin(), exact_counts[0].end(), 0) << endl;
	}
	else {
		apxt_counts = count_approximate(seqs);
		timer.finish();
		verbose_stream << "Seeds = " << apxt_counts[0] << endl;
	}

	timer.go("Allocating hash tables");
	for (unsigned shape_id = shape_from; shape_id < shape_to; ++shape_id)
		assign_ptr(seed_index[shape_id], exact ? new Seed_double_index(exact_counts[shape_id - shape_from]) : new Seed_double_index(apxt_counts[shape_id - shape_from] / Hashed_seed::p));

	timer.go("Building hash table");
	seqs.enum_seeds_partitioned<Count_query_callback, Seed_entry>();

	timer.go("test");
	vector<Access> v(config.threads_);
	seqs.enum_seeds(v);
}
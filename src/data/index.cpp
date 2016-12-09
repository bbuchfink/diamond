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

const double Seed_index::hash_table_factor = 1.3;
Seed_index seed_index[Const::max_shapes];

Seed_index::Seed_index(const Partitioned_histogram &hst, const Sequence_set &seqs, unsigned sid):
	list_buffer(sorted_list::alloc_buffer(hst)),
	list(list_buffer.get(), seqs, shapes[sid], hst.get(sid), seedp_range::all(), hst.partition())
{
	task_timer timer("Counting seeds", 3);
	Thread_pool threads;
	vector<size_t> counts(Const::seedp);
	Atomic<unsigned> seedp(0);
	for (unsigned i = 0; i < config.threads_; ++i)
		threads.push_back(launch_thread(count_seeds, &seedp, &counts, &list));
	threads.join_all();

	timer.go("Allocating hash tables");
	for (unsigned p = 0; p < Hashed_seed::p; ++p)
		assign_ptr(tables[p], new PHash_table<Entry>(counts[p], hash_table_factor));
	
	timer.go("Filling hash tables");
	seedp = 0;
	for (unsigned i = 0; i < config.threads_; ++i)
		threads.push_back(launch_thread(fill_tables, &seedp, this));
	threads.join_all();
}

void Seed_index::count_seeds(Atomic<unsigned> *seedp, vector<size_t> *counts, sorted_list *list)
{
	unsigned p;
	while ((p = (*seedp)++) < Const::seedp) {
		sorted_list::const_iterator i = list->get_partition_cbegin(p);
		size_t n = 0;
		while (!i.at_end()) {
			++n;
			++i;
		}
		(*counts)[p] = n;
	}
}

void Seed_index::fill_tables(Atomic<unsigned> *seedp, Seed_index *idx)
{
	unsigned p;
	while ((p = (*seedp)++) < Const::seedp) {
		sorted_list::const_iterator i = idx->list.get_partition_cbegin(p);
		while (!i.at_end()) {
			PHash_table<Entry>::entry *e = idx->tables[p].insert(murmur_hash()(i.key()));
			e->value.n = (unsigned)idx->list.iterator_offset(i, p);
			++i;
		}
	}
}

void build_index(const Sequence_set &seqs)
{
	task_timer timer("Building histograms", 3);
	const pair<size_t, size_t> len_bounds = seqs.len_bounds(shapes[0].length_);
	const Partitioned_histogram hst(seqs, (unsigned)len_bounds.second);
	timer.finish();

	for (unsigned i = 0; i < shapes.count();++i)
		assign_ptr(seed_index[i], new Seed_index(hst, seqs, i));
}
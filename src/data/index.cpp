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
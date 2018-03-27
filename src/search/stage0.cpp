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

#include "search.h"
#include "../util/algo/hash_join.h"
#include "../util/algo/radix_sort.h"

sorted_list::entry *buf_q;
sorted_list::entry *buf_r;

void seed_join_worker(const sorted_list *query_seeds, const sorted_list *ref_seeds, Atomic<unsigned> *seedp)
{
	unsigned p;
	while ((p = (*seedp)++) < Const::seedp) {
		if (config.sort_join) {
			Relation<sorted_list::entry> R(query_seeds->ptr_begin(p), query_seeds->ptr_end(p) - query_seeds->ptr_begin(p)),
				S(ref_seeds->ptr_begin(p), ref_seeds->ptr_end(p) - ref_seeds->ptr_begin(p));
			//radix_sort(R, 24, &buf_q[query_seeds->limits_[p]]);
			//radix_sort(S, 24, &buf_r[ref_seeds->limits_[p]]);
			radix_sort(R, 24);
			radix_sort(S, 24);
		}
		else {
			hash_join(
				Relation<sorted_list::entry>(query_seeds->ptr_begin(p), query_seeds->ptr_end(p) - query_seeds->ptr_begin(p)),
				Relation<sorted_list::entry>(ref_seeds->ptr_begin(p), ref_seeds->ptr_end(p) - ref_seeds->ptr_begin(p)),
				24
				);
		}
	}
}

void search(const sorted_list &query_seeds, const sorted_list &ref_seeds)
{
	task_timer timer("Allocating radix sort buffer");
	Atomic<unsigned> seedp = 0;
	Thread_pool threads;
	//buf_q = new sorted_list::entry[query_seeds.limits_.back()];
	//buf_r = new sorted_list::entry[ref_seeds.limits_.back()];
	MemoryPool::init(ref_seeds.limits_.back() * sizeof(sorted_list::entry)*config.threads_ * 3 / 1024);
	timer.go("Computing hash join");
	for (size_t i = 0; i < config.threads_; ++i)
		threads.push_back(launch_thread(seed_join_worker, &query_seeds, &ref_seeds, &seedp));
	threads.join_all();
	timer.go("Deallocating radix sort buffer");
	delete[] buf_q;
	delete[] buf_r;
}

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

#include "sorted_list.h"

char* sorted_list::alloc_buffer(const Partitioned_histogram &hst)
{
	return new char[sizeof(entry) * hst.max_chunk_size()];
}

sorted_list::sorted_list()
{}

struct Build_callback
{
	Build_callback(const seedp_range &range, sorted_list::entry* const* ptr, const Seed_set *filter):
		range(range),
		filter(filter),
		it(new sorted_list::buffered_iterator(ptr))
	{ }
	void operator()(uint64_t seed, uint64_t pos, size_t shape)
	{		
		if (filter == 0 || filter->contains(seed))
			it->push(seed, pos, range);
	}
	void finish()
	{
		it->flush();
	}
	seedp_range range;
	const Seed_set *filter;
	auto_ptr<sorted_list::buffered_iterator> it;
};

sorted_list::sorted_list(char *buffer, const Sequence_set &seqs, size_t sh, const shape_histogram &hst, const seedp_range &range, const vector<size_t> seq_partition, const Seed_set *filter) :
	limits_(hst, range),
	data_(reinterpret_cast<entry*>(buffer))
{
	task_timer timer("Building seed list", 3);
	Ptr_set iterators(build_iterators(hst));
	vector<Build_callback> cb(seq_partition.size() - 1);
	for (size_t i = 0; i < seq_partition.size() - 1; ++i)
		cb[i] = Build_callback(range, iterators[i].begin(), filter);

	seqs.enum_seeds(cb, seq_partition, sh, sh + 1);
	timer.go("Sorting seed list");
	Sort_context sort_context(*this);
	launch_scheduled_thread_pool(sort_context, Const::seedp, config.threads_);
}

sorted_list::const_iterator sorted_list::get_partition_cbegin(unsigned p) const
{
	return const_iterator(cptr_begin(p), cptr_end(p));
}

sorted_list::iterator sorted_list::get_partition_begin(unsigned p) const
{
	return iterator(ptr_begin(p), ptr_end(p));
}

sorted_list::Random_access_iterator sorted_list::random_access(unsigned p, size_t offset) const
{
	return Random_access_iterator(cptr_begin(p) + offset, cptr_end(p));
}

sorted_list::entry* sorted_list::ptr_begin(unsigned i) const
{
	return &data_[limits_[i]];
}

sorted_list::entry* sorted_list::ptr_end(unsigned i) const
{
	return &data_[limits_[i + 1]];
}

const sorted_list::entry* sorted_list::cptr_begin(unsigned i) const
{
	return &data_[limits_[i]];
}

const sorted_list::entry* sorted_list::cptr_end(unsigned i) const
{
	return &data_[limits_[i + 1]];
}

void sorted_list::build_seqp(const Sequence_set &seqs, size_t begin, size_t end, entry* const* ptr, const shape &sh, const seedp_range &range, const Seed_set *filter)
{
	uint64_t key;
	auto_ptr<buffered_iterator> it(new buffered_iterator(ptr));
	for (size_t i = begin; i < end; ++i) {
		const sequence seq = seqs[i];
		if (seq.length() < sh.length_) continue;
		const unsigned j1 = (unsigned)seq.length() - sh.length_ + 1;
		for (unsigned j = 0; j < j1; ++j) {
			if (sh.set_seed(key, &seq[j]))
				if (filter == 0 || filter->contains(key))
					it->push(key, seqs.position(i, j), range);
		}
	}
	it->flush();
}

sorted_list::Ptr_set sorted_list::build_iterators(const shape_histogram &hst) const
{
	Ptr_set iterators(hst.size());
	for (unsigned i = 0; i < Const::seedp; ++i)
		iterators[0][i] = ptr_begin(i);

	for (unsigned i = 1; i < hst.size(); ++i)
		for (unsigned j = 0; j < Const::seedp; ++j)
			iterators[i][j] = iterators[i - 1][j] + hst[i - 1][j];
	return iterators;
}
/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

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

#include <stdint.h>
#include "seed_array.h"
#include "seed_set.h"

typedef vector<Array<SeedArray::Entry*, Const::seedp> > PtrSet;

struct BufferedWriter
{
	static const unsigned BUFFER_SIZE = 16;
	BufferedWriter(SeedArray::Entry* const* ptr)
	{
		memset(n, 0, sizeof(n));
		memcpy(this->ptr, ptr, sizeof(this->ptr));
	}
	void push(Packed_seed key, Loc value, const SeedPartitionRange &range)
	{
		const unsigned p = seed_partition(key);
		if (range.contains(p)) {
			assert(n[p] < BUFFER_SIZE);
			buf[p][n[p]++] = SeedArray::Entry(seed_partition_offset(key), value);
			if (n[p] == BUFFER_SIZE)
				flush(p);
		}
	}
	void flush(unsigned p)
	{
		memcpy(ptr[p], buf[p], n[p] * sizeof(SeedArray::Entry));
		ptr[p] += n[p];
		n[p] = 0;
	}
	void flush()
	{
		for (unsigned p = 0; p < Const::seedp; ++p)
			if (n[p] > 0)
				flush(p);
	}
	SeedArray::Entry *ptr[Const::seedp], buf[Const::seedp][BUFFER_SIZE];
	uint8_t n[Const::seedp];
};

PtrSet build_iterators(SeedArray &sa, const shape_histogram &hst)
{
	PtrSet iterators(hst.size());
	for (unsigned i = 0; i < Const::seedp; ++i)
		iterators[0][i] = sa.begin(i);

	for (unsigned i = 1; i < hst.size(); ++i)
		for (unsigned j = 0; j < Const::seedp; ++j)
			iterators[i][j] = iterators[i - 1][j] + hst[i - 1][j];
	return iterators;
}

struct BuildCallback
{
	BuildCallback(const SeedPartitionRange &range, SeedArray::Entry* const* ptr) :
		range(range),
		it(new BufferedWriter(ptr))
	{ }
	bool operator()(uint64_t seed, uint64_t pos, size_t shape)
	{
		it->push(seed, pos, range);
		return true;
	}
	void finish()
	{
		it->flush();
	}
	~BuildCallback()
	{
		delete it;
	}
	SeedPartitionRange range;
	BufferedWriter *it;
};

template<typename _filter>
SeedArray::SeedArray(const Sequence_set &seqs, size_t shape, const shape_histogram &hst, const SeedPartitionRange &range, const vector<size_t> &seq_partition, const _filter *filter):
	range_(range)
{
	task_timer timer("Allocating memory for seed array", 3);
	for (size_t i = range.begin(); i < range.end(); ++i) {
		size_[i] = partition_size(hst, i);
		data_[i] = (Entry*)malloc(size_[i] * sizeof(Entry));
	}
	timer.go("Building seed array");
	PtrSet iterators(build_iterators(*this, hst));
	PtrVector<BuildCallback> cb;
	for (size_t i = 0; i < seq_partition.size() - 1; ++i)
		cb.push_back(new BuildCallback(range, iterators[i].begin()));
	seqs.enum_seeds(cb, seq_partition, shape, shape + 1, filter);
}

SeedArray::~SeedArray() {
	for (size_t i = range_.begin(); i < range_.end(); ++i)
		free(data_[i]);
}

template SeedArray::SeedArray(const Sequence_set &, size_t, const shape_histogram &, const SeedPartitionRange &, const vector<size_t>&, const No_filter *);
template SeedArray::SeedArray(const Sequence_set &, size_t, const shape_histogram &, const SeedPartitionRange &, const vector<size_t>&, const Seed_set *);
template SeedArray::SeedArray(const Sequence_set &, size_t, const shape_histogram &, const SeedPartitionRange &, const vector<size_t>&, const Hashed_seed_set *);
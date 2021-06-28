/****
DIAMOND protein aligner
Copyright (C) 2013-2021 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen
						
Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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
#include "enum_seeds.h"
#include "../util/data_structures/deque.h"

using std::array;

typedef vector<array<SeedArray::Entry*, Const::seedp>> PtrSet;

char* SeedArray::alloc_buffer(const Partitioned_histogram &hst, size_t index_chunks)
{
	return new char[sizeof(Entry) * hst.max_chunk_size(index_chunks)];
}

static const size_t seed_bits(const SeedEncoding code) {
	switch (code) {
	case SeedEncoding::HASHED:
		return sizeof(SeedArray::Entry::Key) * 8;
	case SeedEncoding::SPACED_FACTOR:
		return ceil(shapes[0].weight_ * Reduction::reduction.bit_size_exact()) - Const::seedp_bits;
	case SeedEncoding::CONTIGUOUS:
		return shapes[0].length_ * Reduction::reduction.bit_size() - Const::seedp_bits;
	}
}

struct BufferedWriter
{
	static const unsigned BUFFER_SIZE = 16;
	BufferedWriter(SeedArray::Entry* const* ptr)
	{
		memset(n, 0, sizeof(n));
		memcpy(this->ptr, ptr, sizeof(this->ptr));
	}
	void push(Packed_seed key, Loc value, uint32_t block_id, const SeedPartitionRange &range)
	{
		const unsigned p = seed_partition(key);
		if (range.contains(p)) {
			assert(n[p] < BUFFER_SIZE);
			buf[p][n[p]++] = SeedArray::Entry(seed_partition_offset(key), value, block_id);
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
	bool operator()(uint64_t seed, uint64_t pos, uint32_t block_id, size_t shape)
	{
		it->push(seed, pos, block_id, range);
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
SeedArray::SeedArray(SequenceSet &seqs, size_t shape, const shape_histogram &hst, const SeedPartitionRange &range, const vector<size_t> &seq_partition, char *buffer, const _filter *filter, const SeedEncoding code, const std::vector<bool>* skip) :
	key_bits(seed_bits(code)),
	data_((Entry*)buffer)
{
	begin_[range.begin()] = 0;
	for (size_t i = range.begin(); i < range.end(); ++i)
		begin_[i + 1] = begin_[i] + partition_size(hst, i);

	PtrSet iterators(build_iterators(*this, hst));
	PtrVector<BuildCallback> cb;
	for (size_t i = 0; i < seq_partition.size() - 1; ++i)
		cb.push_back(new BuildCallback(range, iterators[i].data()));
	enum_seeds(&seqs, cb, seq_partition, shape, shape + 1, filter, code, skip);
}

template SeedArray::SeedArray(SequenceSet &, size_t, const shape_histogram &, const SeedPartitionRange &, const vector<size_t>&, char *buffer, const No_filter *, const SeedEncoding, const std::vector<bool>*);
template SeedArray::SeedArray(SequenceSet &, size_t, const shape_histogram &, const SeedPartitionRange &, const vector<size_t>&, char *buffer, const SeedSet *, const SeedEncoding, const std::vector<bool>*);
template SeedArray::SeedArray(SequenceSet &, size_t, const shape_histogram &, const SeedPartitionRange &, const vector<size_t>&, char *buffer, const HashedSeedSet *, const SeedEncoding, const std::vector<bool>*);

struct BufferedWriter2
{
	static const unsigned BUFFER_SIZE = 16;
	BufferedWriter2():
		out(Const::seedp)
	{
		memset(n, 0, sizeof(n));
	}
	void push(Packed_seed key, Loc value, const SeedPartitionRange& range)
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
		out[p].push_back(buf[p], n[p]);
		n[p] = 0;
	}
	void flush()
	{
		for (unsigned p = 0; p < Const::seedp; ++p)
			if (n[p] > 0)
				flush(p);
	}
	vector<Deque<SeedArray::Entry, 15>> out;
	SeedArray::Entry buf[Const::seedp][BUFFER_SIZE];
	uint8_t n[Const::seedp];
};

struct BuildCallback2
{
	BuildCallback2(const SeedPartitionRange& range) :
		range(range),
		it(new BufferedWriter2())
	{ }
	bool operator()(uint64_t seed, uint64_t pos, uint32_t block_id, size_t shape)
	{
		it->push(seed, pos, range);
		return true;
	}
	void finish()
	{
		it->flush();
	}
	~BuildCallback2()
	{
		delete it;
	}
	SeedPartitionRange range;
	BufferedWriter2* it;
};

template<typename _filter>
SeedArray::SeedArray(SequenceSet& seqs, size_t shape, const SeedPartitionRange& range, const _filter* filter, const SeedEncoding code, const std::vector<bool>* skip) :
	key_bits(seed_bits(code)),
	data_(nullptr)
{
	const auto seq_partition = seqs.partition(config.threads_);
	PtrVector<BuildCallback2> cb;
	for (size_t i = 0; i < seq_partition.size() - 1; ++i)
		cb.push_back(new BuildCallback2(range));
	enum_seeds(&seqs, cb, seq_partition, shape, shape + 1, filter, code, skip);

	array<size_t, Const::seedp> counts;
	counts.fill(0);
	for (BuildCallback2* p : cb)
		for (size_t i = 0; i < Const::seedp; ++i)
			counts[i] += p->it->out[i].size();

	for (size_t i = 0; i < Const::seedp; ++i) {
		entries_[i].reserve(counts[i]);
		for (BuildCallback2* p : cb)
			p->it->out[i].move(entries_[i]);
	}
}

template SeedArray::SeedArray(SequenceSet&, size_t, const SeedPartitionRange&, const HashedSeedSet*, const SeedEncoding, const std::vector<bool>*);
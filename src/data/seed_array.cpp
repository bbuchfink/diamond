/****
DIAMOND protein aligner
Copyright (C) 2013-2024 Max Planck Society for the Advancement of Science e.V.
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
#include "../search/seed_complexity.h"

using std::array;
using std::vector;

template<typename SeedLoc> using PtrSet = vector<array<typename SeedArray<SeedLoc>::Entry*, Const::seedp>>;

template<typename SeedLoc>
char* SeedArray<SeedLoc>::alloc_buffer(const SeedHistogram &hst, int index_chunks)
{
	return new char[sizeof(Entry) * hst.max_chunk_size(index_chunks)];
}

template char* SeedArray<PackedLoc>::alloc_buffer(const SeedHistogram&, int);
template char* SeedArray<PackedLocId>::alloc_buffer(const SeedHistogram&, int);

static int seed_bits(const SeedEncoding code) {
	switch (code) {
	case SeedEncoding::HASHED:
		return int(sizeof(SeedOffset) * 8);
	case SeedEncoding::SPACED_FACTOR:
		return int(ceil(shapes[0].weight_ * Reduction::reduction.bit_size_exact()) - Const::seedp_bits);
	case SeedEncoding::CONTIGUOUS:
		return shapes[0].length_ * Reduction::reduction.bit_size() - Const::seedp_bits;
	default:
		break;
	}
	throw std::runtime_error("Unknown seed encoding.");
}

template<typename SeedLoc>
struct BufferedWriter
{
	static const unsigned BUFFER_SIZE = 16;
	BufferedWriter(typename SeedArray<SeedLoc>::Entry* const* ptr)
	{
		memset(n, 0, sizeof(n));
		memcpy(this->ptr, ptr, sizeof(this->ptr));
	}
	void push(PackedSeed key, int64_t value, uint32_t block_id, const SeedPartitionRange &range)
	{
		const unsigned p = seed_partition(key);
		if (range.contains(p)) {
			assert(n[p] < BUFFER_SIZE);
			buf[p][n[p]++] = typename SeedArray<SeedLoc>::Entry(seed_partition_offset(key), value, block_id);
			if (n[p] == BUFFER_SIZE)
				flush(p);
		}
	}
	NO_INLINE void flush(unsigned p)
	{
		memcpy(ptr[p], buf[p], (size_t)n[p] * sizeof(typename SeedArray<SeedLoc>::Entry));
		ptr[p] += n[p];
		n[p] = 0;
	}
	void flush()
	{
		for (unsigned p = 0; p < Const::seedp; ++p)
			if (n[p] > 0)
				flush(p);
	}
	typename SeedArray<SeedLoc>::Entry *ptr[Const::seedp], buf[Const::seedp][BUFFER_SIZE];
	uint8_t n[Const::seedp];
};

template<typename SeedLoc>
PtrSet<SeedLoc> build_iterators(SeedArray<SeedLoc> &sa, const ShapeHistogram &hst)
{
	PtrSet<SeedLoc> iterators(hst.size());
	for (unsigned i = 0; i < Const::seedp; ++i)
		iterators[0][i] = sa.begin(i);

	for (unsigned i = 1; i < hst.size(); ++i)
		for (unsigned j = 0; j < Const::seedp; ++j)
			iterators[i][j] = iterators[i - 1][j] + hst[i - 1][j];
	return iterators;
}

template<typename SeedLoc>
struct BuildCallback
{
	BuildCallback(const SeedPartitionRange &range, typename SeedArray<SeedLoc>::Entry* const* ptr) :
		range(range),
		it(new BufferedWriter<SeedLoc>(ptr))
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
	BufferedWriter<SeedLoc> *it;
};

template<typename SeedLoc> template<typename Filter>
SeedArray<SeedLoc>::SeedArray(Block &seqs, const ShapeHistogram &hst, const SeedPartitionRange &range, char *buffer, const Filter *filter, const EnumCfg& enum_cfg) :
	key_bits(seed_bits(enum_cfg.code)),
	data_((Entry*)buffer)
{
	if (enum_cfg.shape_end - enum_cfg.shape_begin > 1)
		throw std::runtime_error("SeedArray construction for >1 shape.");
	begin_[range.begin()] = 0;
	for (int i = range.begin(); i < range.end(); ++i)
		begin_[i + 1] = begin_[i] + partition_size(hst, i);

	PtrSet<SeedLoc> iterators(build_iterators(*this, hst));
	PtrVector<BuildCallback<SeedLoc>> cb;
	for (size_t i = 0; i < enum_cfg.partition->size() - 1; ++i)
		cb.push_back(new BuildCallback<SeedLoc>(range, iterators[i].data()));
	stats_ = enum_seeds(seqs, cb, filter, enum_cfg);
}

template SeedArray<PackedLoc>::SeedArray(Block&, const ShapeHistogram&, const SeedPartitionRange &, char *buffer, const NoFilter *, const EnumCfg&);
template SeedArray<PackedLoc>::SeedArray(Block&, const ShapeHistogram&, const SeedPartitionRange &, char *buffer, const SeedSet *, const EnumCfg&);
template SeedArray<PackedLoc>::SeedArray(Block&, const ShapeHistogram&, const SeedPartitionRange &, char *buffer, const HashedSeedSet *, const EnumCfg&);
template SeedArray<PackedLocId>::SeedArray(Block&, const ShapeHistogram&, const SeedPartitionRange&, char* buffer, const NoFilter*, const EnumCfg&);
template SeedArray<PackedLocId>::SeedArray(Block&, const ShapeHistogram&, const SeedPartitionRange&, char* buffer, const SeedSet*, const EnumCfg&);
template SeedArray<PackedLocId>::SeedArray(Block&, const ShapeHistogram&, const SeedPartitionRange&, char* buffer, const HashedSeedSet*, const EnumCfg&);

template<typename SeedLoc>
struct BufferedWriter2
{
	static const unsigned BUFFER_SIZE = 16;
	BufferedWriter2():
		out(Const::seedp)
	{
		memset(n, 0, sizeof(n));
	}
	void push(PackedSeed key, int64_t value, const SeedPartitionRange& range)
	{
		const unsigned p = seed_partition(key);
		if (range.contains(p)) {
			assert(n[p] < BUFFER_SIZE);
			buf[p][n[p]++] = typename SeedArray<SeedLoc>::Entry(seed_partition_offset(key), value);
			if (n[p] == BUFFER_SIZE)
				flush(p);
		}
	}
	NO_INLINE void flush(unsigned p)
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
	vector<Deque<typename SeedArray<SeedLoc>::Entry, 15>> out;
	typename SeedArray<SeedLoc>::Entry buf[Const::seedp][BUFFER_SIZE];
	uint8_t n[Const::seedp];
};

template<typename SeedLoc>
struct BuildCallback2
{
	BuildCallback2(const SeedPartitionRange& range) :
		range(range),
		it(new BufferedWriter2<SeedLoc>())
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
	BufferedWriter2<SeedLoc> * it;
};

template<typename SeedLoc> template<typename Filter>
SeedArray<SeedLoc>::SeedArray(Block& seqs, const SeedPartitionRange& range, const Filter* filter, EnumCfg& enum_cfg) :
	key_bits(seed_bits(enum_cfg.code)),
	data_(nullptr)
{
	if (enum_cfg.shape_end - enum_cfg.shape_begin > 1)
		throw std::runtime_error("SeedArray construction for >1 shape.");
	const auto seq_partition = seqs.seqs().partition(config.threads_);
	PtrVector<BuildCallback2<SeedLoc>> cb;
	for (size_t i = 0; i < seq_partition.size() - 1; ++i)
		cb.push_back(new BuildCallback2<SeedLoc>(range));
	enum_cfg.partition = &seq_partition;
	stats_ = enum_seeds(seqs, cb, filter, enum_cfg);

	array<size_t, Const::seedp> counts;
	counts.fill(0);
	for (BuildCallback2<SeedLoc>* p : cb)
		for (size_t i = 0; i < Const::seedp; ++i)
			counts[i] += p->it->out[i].size();

	for (size_t i = 0; i < Const::seedp; ++i) {
		entries_[i].reserve(counts[i]);
		for (BuildCallback2<SeedLoc>* p : cb)
			p->it->out[i].move(entries_[i]);
	}
}

template SeedArray<PackedLoc>::SeedArray(Block&, const SeedPartitionRange&, const HashedSeedSet*, EnumCfg&);
template SeedArray<PackedLocId>::SeedArray(Block&, const SeedPartitionRange&, const HashedSeedSet*, EnumCfg&);
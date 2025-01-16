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
#include "util/data_structures/deque.h"
#include "util/memory/alignment.h"

using std::vector;
using std::fill;
using std::array;

template<typename SeedLoc> using PtrSet = vector<vector<typename SeedArray<SeedLoc>::Entry*>>;

template<typename SeedLoc>
char* SeedArray<SeedLoc>::alloc_buffer(const SeedHistogram &hst, int index_chunks)
{
	return (char*)Util::Memory::aligned_malloc(sizeof(Entry) * hst.max_chunk_size(index_chunks), 32);
	//return new char[sizeof(Entry) * hst.max_chunk_size(index_chunks)];
}

template char* SeedArray<PackedLoc>::alloc_buffer(const SeedHistogram&, int);
template char* SeedArray<PackedLocId>::alloc_buffer(const SeedHistogram&, int);

static int seed_bits(const SeedEncoding code, int seedp_bits) {
	switch (code) {
	case SeedEncoding::HASHED:
		return int(sizeof(SeedOffset) * 8);
	case SeedEncoding::SPACED_FACTOR:
		return int(ceil(shapes[0].weight_ * Reduction::reduction.bit_size_exact()) - seedp_bits);
	case SeedEncoding::CONTIGUOUS:
		return shapes[0].length_ * Reduction::reduction.bit_size() - seedp_bits;
	default:
		break;
	}
	throw std::runtime_error("Unknown seed encoding.");
}

template<typename SeedLoc>
struct BufferedWriter
{
	using Entry = typename SeedArray<SeedLoc>::Entry;
	static const unsigned BUFFER_SIZE = 16;
	BufferedWriter(Entry* const* ptr, int seedp_bits, SeedPartitionRange seedp_range):
		seedp_mask(::seedp_mask(seedp_bits)),
		seedp_bits(seedp_bits),
		ptr(ptr, ptr + seedp_range.size()),
		buf(seedp_range.size()),
		n(seedp_range.size(), 0)
	{}
	void push(PackedSeed key, int64_t value, uint32_t block_id, const SeedPartitionRange &range)
	{
		const SeedPartition p = seed_partition(key, seedp_mask);
		if (range.contains(p)) {
			const SeedPartition d = p - range.begin();
			assert(n[d] < BUFFER_SIZE);
			buf[d][n[d]++] = Entry(seed_partition_offset(key, seedp_bits), value, block_id);
			if (n[d] == BUFFER_SIZE)
				flush(d);
		}
	}
	NO_INLINE void flush(SeedPartition p)
	{
		memcpy(ptr[p], buf[p].data(), (size_t)n[p] * sizeof(Entry));
		ptr[p] += n[p];
		n[p] = 0;
	}
	void flush()
	{
		for (SeedPartition p = 0; p < (SeedPartition)n.size(); ++p)
			if (n[p] > 0)
				flush(p);
	}
	const PackedSeed seedp_mask, seedp_bits;
	vector<Entry*> ptr;
	vector<array<Entry, BUFFER_SIZE>> buf;
	vector<uint8_t> n;
};

template<typename SeedLoc>
static PtrSet<SeedLoc> build_iterators(SeedArray<SeedLoc> &sa, const ShapeHistogram &hst, SeedPartitionRange range)
{
	PtrSet<SeedLoc> iterators(hst.size());
	for (auto& v : iterators)
		v.reserve(range.size());
	for (SeedPartition i = 0; i < range.size(); ++i)
		iterators[0].push_back(sa.begin(i));

	for (unsigned i = 1; i < hst.size(); ++i)
		for (SeedPartition j = range.begin(); j < range.end(); ++j)
			iterators[i].push_back(iterators[i - 1][j - range.begin()] + hst[i - 1][j]);
	return iterators;
}

template<typename SeedLoc>
struct BuildCallback
{
	BuildCallback(const SeedPartitionRange &range, typename SeedArray<SeedLoc>::Entry* const* ptr, int seedp_bits) :
		range(range),
		it(new BufferedWriter<SeedLoc>(ptr, seedp_bits, range))
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
	const SeedPartitionRange range;
	BufferedWriter<SeedLoc> *it;
};

template<typename SeedLoc> template<typename Filter>
SeedArray<SeedLoc>::SeedArray(Block &seqs, const ShapeHistogram &hst, const SeedPartitionRange &range, int seedp_bits, char *buffer, const Filter *filter, const EnumCfg& enum_cfg) :
	key_bits(seed_bits(enum_cfg.code, seedp_bits)),
	data_((Entry*)buffer)
{
	if (enum_cfg.shape_end - enum_cfg.shape_begin > 1)
		throw std::runtime_error("SeedArray construction for >1 shape.");
	begin_.reserve(range.size() + 1);
	begin_.push_back(0);
	for (SeedPartition i = range.begin(); i < range.end(); ++i)
		begin_.push_back(begin_.back() + partition_size(hst, i));

	PtrSet<SeedLoc> iterators(build_iterators(*this, hst, range));
	PtrVector<BuildCallback<SeedLoc>> cb;
	for (size_t i = 0; i < enum_cfg.partition->size() - 1; ++i)
		cb.push_back(new BuildCallback<SeedLoc>(range, iterators[i].data(), seedp_bits));
	stats_ = enum_seeds(seqs, cb, filter, enum_cfg);
}

template SeedArray<PackedLoc>::SeedArray(Block&, const ShapeHistogram&, const SeedPartitionRange&, int, char* buffer, const NoFilter*, const EnumCfg&);
template SeedArray<PackedLoc>::SeedArray(Block&, const ShapeHistogram&, const SeedPartitionRange&, int, char* buffer, const SeedSet*, const EnumCfg&);
template SeedArray<PackedLoc>::SeedArray(Block&, const ShapeHistogram&, const SeedPartitionRange&, int, char* buffer, const HashedSeedSet*, const EnumCfg&);
template SeedArray<PackedLocId>::SeedArray(Block&, const ShapeHistogram&, const SeedPartitionRange&, int, char* buffer, const NoFilter*, const EnumCfg&);
template SeedArray<PackedLocId>::SeedArray(Block&, const ShapeHistogram&, const SeedPartitionRange&, int, char* buffer, const SeedSet*, const EnumCfg&);
template SeedArray<PackedLocId>::SeedArray(Block&, const ShapeHistogram&, const SeedPartitionRange&, int, char* buffer, const HashedSeedSet*, const EnumCfg&);

template<typename SeedLoc>
struct OnePassBufferedWriter
{
	using Entry = typename SeedArray<SeedLoc>::Entry;
	static const unsigned BUFFER_SIZE = 16;
	OnePassBufferedWriter(SeedPartitionRange range, int seedp_bits) :
		seedp_bits(seedp_bits),
		out(range.size()),
		buf(range.size()),
		n(range.size(), 0)
	{}
	void push(PackedSeed key, int64_t value, const SeedPartitionRange& range)
	{
		const SeedPartition p = seed_partition(key, seedp_mask(seedp_bits));
		if (range.contains(p)) {
			const SeedPartition d = p - range.begin();
			assert(n[d] < BUFFER_SIZE);
			buf[d][n[d]++] = Entry(seed_partition_offset(key, seedp_bits), value);
			if (n[d] == BUFFER_SIZE)
				flush(d);
		}
	}
	NO_INLINE void flush(unsigned p)
	{
		out[p].push_back(buf[p].data(), n[p]);
		n[p] = 0;
	}
	void flush()
	{
		for (SeedPartition p = 0; p < (SeedPartition)n.size(); ++p)
			if (n[p] > 0)
				flush(p);
	}
	const int seedp_bits;
	vector<Deque<Entry, 15>> out;
	vector<array<Entry, BUFFER_SIZE>> buf;
	vector<uint8_t> n;
};

template<typename SeedLoc>
struct OnePassBuildCallback
{
	OnePassBuildCallback(const SeedPartitionRange& range, int seedp_bits) :
		range(range),
		it(new OnePassBufferedWriter<SeedLoc>(range, seedp_bits))
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
	~OnePassBuildCallback()
	{
		delete it;
	}
	SeedPartitionRange range;
	OnePassBufferedWriter<SeedLoc>* it;
};

template<typename SeedLoc> template<typename Filter>
SeedArray<SeedLoc>::SeedArray(Block& seqs, const SeedPartitionRange& range, int seedp_bits, const Filter* filter, EnumCfg& enum_cfg) :
	key_bits(seed_bits(enum_cfg.code, seedp_bits)),
	data_(nullptr)
{
	if (enum_cfg.shape_end - enum_cfg.shape_begin > 1)
		throw std::runtime_error("SeedArray construction for >1 shape.");
	const auto seq_partition = seqs.seqs().partition(config.threads_);
	PtrVector<OnePassBuildCallback<SeedLoc>> cb;
	for (size_t i = 0; i < seq_partition.size() - 1; ++i)
		cb.push_back(new OnePassBuildCallback<SeedLoc>(range, seedp_bits));
	enum_cfg.partition = &seq_partition;
	stats_ = enum_seeds(seqs, cb, filter, enum_cfg);

	vector<size_t> counts(range.size(), 0);
	for (OnePassBuildCallback<SeedLoc>* p : cb)
		for (SeedPartition i = 0; i < range.size(); ++i)
			counts[i] += p->it->out[i].size();

	for (SeedPartition i = 0; i < range.size(); ++i) {
		entries_[i].reserve(counts[i]);
		for (OnePassBuildCallback<SeedLoc>* p : cb)
			p->it->out[i].move(entries_[i]);
	}
}

template SeedArray<PackedLoc>::SeedArray(Block&, const SeedPartitionRange&, int, const HashedSeedSet*, EnumCfg&);
template SeedArray<PackedLocId>::SeedArray(Block&, const SeedPartitionRange&, int, const HashedSeedSet*, EnumCfg&);
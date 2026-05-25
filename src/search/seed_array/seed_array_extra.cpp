/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

#include <stdint.h>
#include "seed_array.h"
#include "data/seed_set.h"
#include "enum_seeds.h"
#include "util/data_structures/deque.h"
#include "seed_array_impl.h"

using std::vector;
using std::fill;
using std::array;

namespace DISPATCH_ARCH {

template char* SeedArray<PackedLocId>::alloc_buffer(const SeedHistogram&, int);

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
		seedp_mask(::seedp_mask(seedp_bits)),
		seedp_bits(seedp_bits),
		range_begin(range.begin()),
		range_size(range.size()),
		out(range.size()),
		buf(range.size()),
		n(range.size(), 0)
	{
	}
	void push(PackedSeed key, int64_t value)
	{
		const SeedPartition d = seed_partition(key, seedp_mask) - range_begin;
		if (d < range_size) {
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
	const PackedSeed seedp_mask;
	const int seedp_bits;
	const SeedPartition range_begin, range_size;
	vector<Deque<Entry, 15>> out;
	vector<array<Entry, BUFFER_SIZE>> buf;
	vector<uint8_t> n;
};

template<typename SeedLoc>
struct OnePassBuildCallback
{
	OnePassBuildCallback(const SeedPartitionRange& range, int seedp_bits) :
		it(range, seedp_bits)
	{
	}
	bool operator()(uint64_t seed, uint64_t pos, uint32_t block_id, size_t shape)
	{
		it.push(seed, pos);
		return true;
	}
	void finish()
	{
		it.flush();
	}
	OnePassBufferedWriter<SeedLoc> it;
};

template<typename SeedLoc> template<typename Filter>
SeedArray<SeedLoc>::SeedArray(Block& seqs, const SeedPartitionRange& range, int seedp_bits, const Filter* filter, EnumCfg& enum_cfg) :
	key_bits(seed_bits(enum_cfg.code, seedp_bits)),
	data_(nullptr),
	entries_(range.size())
{
	if (enum_cfg.shape_end - enum_cfg.shape_begin > 1)
		throw std::runtime_error("SeedArray construction for >1 shape.");
	const auto seq_partition = seqs.seqs().partition(config.threads_);
	PtrVector<OnePassBuildCallback<SeedLoc>> cb;
	cb.reserve(seq_partition.size() - 1);
	for (size_t i = 0; i < seq_partition.size() - 1; ++i)
		cb.push_back(new OnePassBuildCallback<SeedLoc>(range, seedp_bits));
	enum_cfg.partition = &seq_partition;
	stats_ = enum_seeds(seqs, cb, filter, enum_cfg);

	vector<size_t> counts(range.size(), 0);
	for (OnePassBuildCallback<SeedLoc>* p : cb)
		for (SeedPartition i = 0; i < range.size(); ++i)
			counts[i] += p->it.out[i].size();

	for (SeedPartition i = 0; i < range.size(); ++i) {
		entries_[i].reserve(counts[i]);
		for (OnePassBuildCallback<SeedLoc>* p : cb)
			p->it.out[i].move(entries_[i]);
	}
}

template SeedArray<PackedLoc>::SeedArray(Block&, const SeedPartitionRange&, int, const HashedSeedSet*, EnumCfg&);
template SeedArray<PackedLocId>::SeedArray(Block&, const SeedPartitionRange&, int, const HashedSeedSet*, EnumCfg&);

}

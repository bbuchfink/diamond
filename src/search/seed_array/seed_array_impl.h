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

#pragma once
#include "seed_array.h"

using std::vector;
using std::array;

namespace DISPATCH_ARCH {

static inline int seed_bits(const SeedEncoding code, int seedp_bits) {
	switch (code) {
	case SeedEncoding::HASHED:
		return int(sizeof(SeedOffset) * 8);
	case SeedEncoding::SPACED_FACTOR:
		return int(ceil(shapes[0].weight_ * Reduction::get_reduction().bit_size_exact()) - seedp_bits);
	case SeedEncoding::CONTIGUOUS:
		return shapes[0].length_ * Reduction::get_reduction().bit_size() - seedp_bits;
	default:
		break;
	}
	throw std::runtime_error("Unknown seed encoding.");
}

template<typename SeedLoc>
struct BufferedWriter
{
	using Entry = typename SeedArray<SeedLoc>::Entry;
	static const unsigned BUFFER_SIZE = 64;
	struct PartitionState {
		uint32_t n;
		Entry buf[BUFFER_SIZE];
		Entry* ptr;
	};
	BufferedWriter(Entry* const* ptr, int seedp_bits, SeedPartitionRange seedp_range) :
		seedp_mask(::seedp_mask(seedp_bits)),
		seedp_bits(seedp_bits),
		range_begin(seedp_range.begin()),
		range_size(seedp_range.size()),
		state(seedp_range.size())
	{
		for (SeedPartition i = 0; i < range_size; ++i) {
			state[i].n = 0;
			state[i].ptr = ptr[i];
		}
	}
	void push(PackedSeed key, int64_t value, uint32_t block_id)
	{
		const SeedPartition p = seed_partition(key, seedp_mask) - range_begin;
		if (p < range_size) {
			PartitionState& s = state[p];
			assert(s.n < BUFFER_SIZE);
			s.buf[s.n++] = Entry(seed_partition_offset(key, seedp_bits), value, block_id);
			if (s.n == BUFFER_SIZE)
				flush(s);
		}
	}
	NO_INLINE void flush(PartitionState& s)
	{
		memcpy(s.ptr, s.buf, (size_t)s.n * sizeof(Entry));
		s.ptr += s.n;
		s.n = 0;
	}
	void flush()
	{
		for (auto& s : state)
			if (s.n > 0)
				flush(s);
	}
	const PackedSeed seedp_mask;
	const int seedp_bits;
	const SeedPartition range_begin, range_size;
	vector<PartitionState> state;
};

template<typename SeedLoc>
struct BuildCallback
{
	BuildCallback(const SeedPartitionRange& range, typename SeedArray<SeedLoc>::Entry* const* ptr, int seedp_bits) :
		it(ptr, seedp_bits, range)
	{
	}
	bool operator()(uint64_t seed, uint64_t pos, uint32_t block_id, size_t shape)
	{
		it.push(seed, pos, block_id);
		return true;
	}
	void finish()
	{
		it.flush();
	}
	BufferedWriter<SeedLoc> it;
};

template<typename SeedLoc> using PtrSet = std::vector<std::vector<typename SeedArray<SeedLoc>::Entry*>>;

template<typename SeedLoc>
static PtrSet<SeedLoc> build_iterators(SeedArray<SeedLoc>& sa, const ShapeHistogram& hst, SeedPartitionRange range)
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

template<typename SeedLoc> template<typename Filter>
SeedArray<SeedLoc>::SeedArray(Block& seqs, const ShapeHistogram& hst, const SeedPartitionRange& range, int seedp_bits, char* buffer, const Filter* filter, const EnumCfg& enum_cfg) :
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
	cb.reserve(enum_cfg.partition->size() - 1);
	for (size_t i = 0; i < enum_cfg.partition->size() - 1; ++i)
		cb.push_back(new BuildCallback<SeedLoc>(range, iterators[i].data(), seedp_bits));
	stats_ = enum_seeds(seqs, cb, filter, enum_cfg);
}

}
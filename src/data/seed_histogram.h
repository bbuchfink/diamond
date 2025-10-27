/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <vector>
#include "flags.h"
#include "basic/seed.h"

using ShapeHistogram = std::vector<std::vector<unsigned>>;

struct SeedPartitionRange
{
	SeedPartitionRange():
		begin_ (0),
		end_ (0)
	{ }
	SeedPartitionRange(SeedPartition begin, SeedPartition end):
		begin_ (begin),
		end_ (end)
	{ }
	bool contains(SeedPartition i) const
	{ return i >= begin_ && i < end_; }
	SeedPartition begin() const
	{ return begin_; }
	SeedPartition end() const
	{ return end_; }
	bool lower(SeedPartition i) const
	{ return i < begin_; }
	bool lower_or_equal(SeedPartition i) const
	{ return i < end_; }
	SeedPartition size() const
	{
		return end_ - begin_;
	}
private:
	SeedPartition begin_, end_;
};

extern SeedPartitionRange current_range;

inline size_t partition_size(const ShapeHistogram &hst, size_t p)
{
	size_t s = 0;
	for (unsigned i = 0; i < hst.size(); ++i)
		s += hst[i][p];
	return s;
}

inline size_t hst_size(const ShapeHistogram &hst, const SeedPartitionRange &range)
{
	size_t s = 0;
	for (SeedPartition i = range.begin(); i < range.end(); ++i)
		s += partition_size(hst, i);
	return s;
}

struct Block;

struct SeedHistogram
{

	SeedHistogram();
	
	template<typename Filter>
	SeedHistogram(Block& seqs, bool serial, const Filter* filter, EnumCfg& enum_cfg, int seedp_bits);

	const ShapeHistogram& get(unsigned sid) const
	{ return data_[sid]; }

	size_t max_chunk_size(const int index_chunks) const;

	const std::vector<uint32_t>& partition() const
	{
		return p_;
	}

	int seedp() const {
		return (int)data_.front().front().size();
	}

private:

	std::vector<uint32_t> p_;
	std::vector<ShapeHistogram> data_;

};
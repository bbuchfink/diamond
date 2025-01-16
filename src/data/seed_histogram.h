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

	const std::vector<BlockId>& partition() const
	{
		return p_;
	}

	int seedp() const {
		return (int)data_.front().front().size();
	}

private:

	std::vector<BlockId> p_;
	std::vector<ShapeHistogram> data_;

};
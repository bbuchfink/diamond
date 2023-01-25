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
#include <array>
#include <vector>
#include "flags.h"
#include "../basic/const.h"
#include "../masking/masking.h"

typedef std::vector<std::array<unsigned, Const::seedp>> ShapeHistogram;

struct SeedPartitionRange
{
	SeedPartitionRange():
		begin_ (0),
		end_ (0)
	{ }
	SeedPartitionRange(int begin, int end):
		begin_ (begin),
		end_ (end)
	{ }
	bool contains(int i) const
	{ return i >= begin_ && i < end_; }
	int begin() const
	{ return begin_; }
	int end() const
	{ return end_; }
	bool lower(int i) const
	{ return i < begin_; }
	bool lower_or_equal(int i) const
	{ return i < end_; }
	int size() const
	{
		return end_ - begin_;
	}
	static SeedPartitionRange all()
	{
		return SeedPartitionRange(0, Const::seedp);
	}
private:
	int begin_, end_;
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
	for(int i=range.begin();i<range.end();++i)
		s += partition_size(hst, i);
	return s;
}

struct Block;

struct SeedHistogram
{

	SeedHistogram();
	
	template<typename Filter>
	SeedHistogram(Block& seqs, bool serial, const Filter* filter, SeedEncoding code, const std::vector<bool>* skip, const bool mask_seeds, const double seed_cut, const MaskingAlgo soft_masking, Loc minimizer_window);

	const ShapeHistogram& get(unsigned sid) const
	{ return data_[sid]; }

	size_t max_chunk_size(const int index_chunks) const;

	const std::vector<size_t>& partition() const
	{
		return p_;
	}

private:

	std::vector<ShapeHistogram> data_;
	std::vector<size_t> p_;

};
/****
DIAMOND protein aligner
Copyright (C) 2016-2021 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
						
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
#include "../util/ptr_vector.h"
#include "../util/data_structures/hash_set.h"
#include "../masking/masking.h"
#include "../lib/mio/forward.h"

const uint64_t SEED_INDEX_MAGIC_NUMBER = 0x2d6ba306ecbf6aba;
const uint32_t SEED_INDEX_VERSION = 0;
const size_t SEED_INDEX_HEADER_SIZE = 16;

struct Block;

struct SeedSet
{
	SeedSet(Block &seqs, double max_coverage, const std::vector<bool>* skip, const double seed_cut, const MaskingAlgo soft_masking);
	bool contains(uint64_t key, uint64_t shape) const
	{
		return data_[key];
	}
	double coverage() const
	{
		return coverage_;
	}
private:
	std::vector<bool> data_;
	double coverage_;
};

struct HashedSeedSet
{
	typedef HashSet<Modulo2, Identity> Table;
	HashedSeedSet(Block &seqs, const std::vector<bool>* skip, const double seed_cut, const MaskingAlgo soft_masking);
	HashedSeedSet(const std::string& index_file);
	~HashedSeedSet();
	bool contains(uint64_t key, uint64_t shape) const
	{
		return data_[shape].contains(key);
	}
	const Table& table(size_t i) const {
		return data_[i];
	}
	size_t max_table_size() const;
private:
	PtrVector<Table> data_;
	std::unique_ptr<mio::mmap_source> mmap_;
};

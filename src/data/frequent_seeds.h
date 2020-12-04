/****
DIAMOND protein aligner
Copyright (C) 2016-2020 Max Planck Society for the Advancement of Science e.V.
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
#include <atomic>
#include "../basic/const.h"
#include "../util/hash_table.h"
#include "seed_array.h"
#include "../util/algo/join_result.h"
#include "../util/range.h"

struct Frequent_seeds
{

	void build(unsigned sid, const SeedPartitionRange &range, DoubleArray<SeedArray::_pos> *query_seed_hits, DoubleArray<SeedArray::_pos> *ref_seed_hits);

	bool get(const Letter *pos, unsigned sid) const
	{
		Packed_seed seed;
		const bool t = config.algo == Config::double_indexed ? shapes[sid].set_seed(seed, pos) : shapes[sid].set_seed_shifted(seed, pos);
		if (!t)
			return true;
		return tables_[sid][seed_partition(seed)].contains(seed_partition_offset(seed));
	}

	static void clear_masking(Sequence_set& seqs);

private:

	static const double hash_table_factor;   

	static void build_worker(
		size_t seedp,
		size_t thread_id,
		DoubleArray<SeedArray::_pos> *query_seed_hits,
		DoubleArray<SeedArray::_pos> *ref_seed_hits,
		const SeedPartitionRange *range,
		unsigned sid,
		unsigned ref_max_n,
		unsigned query_max_n,
		vector<unsigned> *counts);

	static void compute_sd(std::atomic<unsigned> *seedp, DoubleArray<SeedArray::_pos> *query_seed_hits, DoubleArray<SeedArray::_pos> *ref_seed_hits, vector<Sd> *ref_out, vector<Sd> *query_out);

	PHash_set<void,murmur_hash> tables_[Const::max_shapes][Const::seedp];

};

extern Frequent_seeds frequent_seeds;


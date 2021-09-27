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
#include "flags.h"
#include "../util/algo/join_result.h"
#include "../util/range.h"
#include "../run/config.h"
#include "../util/data_structures/double_array.h"
#include "seed_histogram.h"

struct SequenceSet;

struct Frequent_seeds
{

	void build(unsigned sid, const SeedPartitionRange &range, DoubleArray<SeedLoc> *query_seed_hits, DoubleArray<SeedLoc> *ref_seed_hits, Search::Config& cfg);
	static void clear_masking(SequenceSet& seqs);

private:

	static const double hash_table_factor;   

	static void build_worker(
		size_t seedp,
		size_t thread_id,
		DoubleArray<SeedLoc> *query_seed_hits,
		DoubleArray<SeedLoc> *ref_seed_hits,
		const SeedPartitionRange *range,
		unsigned sid,
		unsigned ref_max_n,
		unsigned query_max_n,
		std::vector<unsigned> *counts,
		Search::Config* cfg);

};

extern Frequent_seeds frequent_seeds;
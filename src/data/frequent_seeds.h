/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#ifndef FREQUENT_SEEDS_H_
#define FREQUENT_SEEDS_H_

#include "../basic/const.h"
#include "../util/hash_table.h"
#include "sorted_list.h"

struct Frequent_seeds
{

	void build(unsigned sid, const seedp_range &range, sorted_list &ref_idx, const sorted_list &query_idx);

	bool get(const Letter *pos, unsigned sid) const
	{
		Packed_seed seed;
		const bool t = config.algo == Config::double_indexed ? shapes[sid].set_seed(seed, pos) : shapes[sid].set_seed_shifted(seed, pos);
		if (!t)
			return true;
		return tables_[sid][seed_partition(seed)].contains(seed_partition_offset(seed));
	}

private:

	static const double hash_table_factor;   

	struct Build_context;

	static void compute_sd(Atomic<unsigned> *seedp, const sorted_list *ref_idx, const sorted_list *query_idx, vector<Sd> *ref_out, vector<Sd> *query_out);

	PHash_set<void,murmur_hash> tables_[Const::max_shapes][Const::seedp];

};

extern Frequent_seeds frequent_seeds;

#endif
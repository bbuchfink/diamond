/****
Copyright (c) 2016, Benjamin Buchfink
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
		if (!shapes[sid].set_seed(seed, pos))
			return true;
		return tables_[sid][seed_partition(seed)].contains(seed_partition_offset(seed));
	}

private:

	static const double hash_table_factor;   

	struct Build_context;

	static void compute_sd(Atomic<unsigned> *seedp, const sorted_list *ref_idx, const sorted_list *query_idx, vector<Sd> *ref_out, vector<Sd> *query_out);

	PHash_set tables_[Const::max_shapes][Const::seedp];

};

extern Frequent_seeds frequent_seeds;

#endif
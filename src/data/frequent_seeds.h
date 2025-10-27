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
#include "run/config.h"
#include "util/data_structures/double_array.h"
#include "seed_histogram.h"

struct SequenceSet;

struct FrequentSeeds
{

	template<typename SeedLoc>
	void build(unsigned sid, const SeedPartitionRange &range, DoubleArray<SeedLoc> *query_seed_hits, DoubleArray<SeedLoc> *ref_seed_hits, Search::Config& cfg);
	static void clear_masking(SequenceSet& seqs);

private:

	static const double hash_table_factor;   

	template<typename SeedLoc>
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

extern FrequentSeeds frequent_seeds;
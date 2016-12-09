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

#include "../util/thread.h"
#include "../data/queries.h"
#include "../data/index.h"

void search_query(unsigned query_id, Statistics &stat, vector<Seed> &neighbor_seeds)
{
	const sequence query_seq = query_seqs::get()[query_id];
	for (unsigned sid = 0; sid < shapes.count(); ++sid) {
		const shape &sh = shapes[sid];
		if (query_seq.length() < sh.length_)
			return;
		for (unsigned i = 0; i <= query_seq.length() - sh.length_; ++i) {
			uint64_t seed;
			if (sh.set_seed(seed, &query_seq[i])) {
				sorted_list::Random_access_iterator k = seed_index[sid][seed];
				while (k.good()) {
					stat.inc(Statistics::SEED_HITS);
					++k;
				}
			}
		}
	}
}

void search_query_worker(Atomic<unsigned> *next)
{
	unsigned query_id;
	Statistics stat;
	vector<Seed> neighbor_seeds;
	while ((query_id = (*next)++) < query_seqs::get().get_length())
		search_query(query_id, stat, neighbor_seeds);
	statistics += stat;
}
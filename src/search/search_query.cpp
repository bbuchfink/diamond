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
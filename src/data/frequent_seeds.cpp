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

#include <numeric>
#include "frequent_seeds.h"
#include "sorted_list.h"

const double Frequent_seeds::hash_table_factor = 1.3;
Frequent_seeds frequent_seeds;

void Frequent_seeds::compute_sd(Atomic<unsigned> *seedp, const sorted_list *ref_idx, const sorted_list *query_idx, vector<Sd> *ref_out, vector<Sd> *query_out)
{
	unsigned p;
	while ((p = (*seedp)++) < current_range.end()) {
		Sd ref_sd;
		sorted_list::const_iterator it = ref_idx->get_partition_cbegin(p);
		while (!it.at_end()) {
			ref_sd.add((double)it.n);
			++it;
		}
		(*ref_out)[p - current_range.begin()] = ref_sd;

		Sd query_sd;
		it = query_idx->get_partition_cbegin(p);
		while (!it.at_end()) {
			query_sd.add((double)it.n);
			++it;
		}
		(*query_out)[p - current_range.begin()] = query_sd;
	}
}

struct Frequent_seeds::Build_context
{
	Build_context(const sorted_list &ref_idx, const sorted_list &query_idx, const seedp_range &range, unsigned sid, unsigned ref_max_n, unsigned query_max_n, vector<unsigned> &counts) :
		ref_idx(ref_idx),
		query_idx(query_idx),
		range(range),
		sid(sid),
		ref_max_n(ref_max_n),
		query_max_n(query_max_n),
		counts(counts)
	{ }
	void operator()(unsigned thread_id, unsigned seedp)
	{
		if (!range.contains(seedp))
			return;
		
		vector<uint32_t> buf;
		size_t n = 0;
		Merge_iterator<sorted_list::iterator> merge_it(ref_idx.get_partition_begin(seedp), query_idx.get_partition_begin(seedp));
		while (merge_it.next()) {
			if (merge_it.i.n > ref_max_n || merge_it.j.n > query_max_n) {
				merge_it.i.get(0)->value = 0;
				n += (unsigned)merge_it.i.n;
				buf.push_back(merge_it.i.key());
			}
			++merge_it;
		}

		const size_t ht_size = std::max((size_t)(buf.size() * hash_table_factor), buf.size() + 1);
		PHash_set<void, murmur_hash> hash_set(ht_size);

		for (vector<uint32_t>::const_iterator i = buf.begin(); i != buf.end(); ++i)
			hash_set.insert(*i);

		frequent_seeds.tables_[sid][seedp] = hash_set;
		counts[seedp] = (unsigned)n;
	}
	const sorted_list &ref_idx;
	const sorted_list &query_idx;
	const seedp_range range;
	const unsigned sid, ref_max_n, query_max_n;
	vector<unsigned> &counts;
};

void Frequent_seeds::build(unsigned sid, const seedp_range &range, sorted_list &ref_idx, const sorted_list &query_idx)
{
	vector<Sd> ref_sds(range.size()), query_sds(range.size());
	Atomic<unsigned> seedp (range.begin());
	Thread_pool threads;
	for (unsigned i = 0; i < config.threads_; ++i)
		threads.push_back(launch_thread(compute_sd, &seedp, &ref_idx, &query_idx, &ref_sds, &query_sds));
	threads.join_all();

	Sd ref_sd(ref_sds), query_sd(query_sds);
	const unsigned ref_max_n = (unsigned)(ref_sd.mean() + config.freq_sd*ref_sd.sd()), query_max_n = (unsigned)(query_sd.mean() + config.freq_sd*query_sd.sd());
	log_stream << "Seed frequency mean (reference) = " << ref_sd.mean() << ", SD = " << ref_sd.sd() << endl;
	log_stream << "Seed frequency mean (query) = " << query_sd.mean() << ", SD = " << query_sd.sd() << endl;
	vector<unsigned> counts(Const::seedp);
	Build_context build_context(ref_idx, query_idx, range, sid, ref_max_n, query_max_n, counts);
	launch_scheduled_thread_pool(build_context, Const::seedp, config.threads_);
	log_stream << "Masked positions = " << std::accumulate(counts.begin(), counts.end(), 0) << std::endl;
}
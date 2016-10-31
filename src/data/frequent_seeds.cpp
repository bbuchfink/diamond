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
		PHash_set hash_set(ht_size);

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
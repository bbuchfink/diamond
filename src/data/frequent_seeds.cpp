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

struct Frequent_seeds::Build_context
{
	Build_context(const sorted_list &ref_idx, const sorted_list &query_idx, unsigned sid, vector<unsigned> &counts) :
		ref_idx(ref_idx),
		query_idx(query_idx),
		sid(sid),
		counts(counts)
	{ }
	void operator()(unsigned thread_id, unsigned seedp)
	{
		Sd ref_sd, query_sd, mult_sd;
		Merge_iterator<sorted_list::const_iterator> it(ref_idx.get_partition_cbegin(seedp), query_idx.get_partition_cbegin(seedp));
		while(it.next()) {
			ref_sd.add((double)it.i.n);
			query_sd.add((double)it.j.n);
			mult_sd.add((double)it.i.n*(double)it.j.n);
			++it;
		}

		/*sorted_list::const_iterator it = ref_idx.get_partition_cbegin(seedp);
		while (!it.at_end()) {
			ref_sd.add((double)it.n);
			++it;
		}

		const unsigned max_n = (unsigned)(ref_sd.mean() + config.freq_sd*ref_sd.sd());
		//cout << "mean=" << A << " sd=" << sd << endl;

		/*size_t n = 0;
		vector<uint32_t> buf;
		sorted_list::iterator j = ref_idx.get_partition_begin(seedp);
		while (!j.at_end()) {
			if (j.n > max_n) {
				j.get(0)->value = 0;
				n += (unsigned)j.n;
				buf.push_back(j.key());
			}
			++j;
		}*/

		const size_t max_n = (unsigned)(mult_sd.mean() + config.freq_sd*mult_sd.sd());
		vector<uint32_t> buf;
		size_t n = 0;
		{
			Merge_iterator<sorted_list::iterator> it(ref_idx.get_partition_begin(seedp), query_idx.get_partition_begin(seedp));
			while (it.next()) {
				if ((size_t)it.i.n * (size_t)it.j.n > max_n) {
				//if (it.i.n > max_n) {
					it.i.get(0)->value = 0;
					n += (unsigned)it.i.n;
					buf.push_back(it.i.key());
				}
				++it;
			}
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
	const unsigned sid;
	vector<unsigned> &counts;
};

void Frequent_seeds::build(unsigned sid, const seedp_range &range, sorted_list &ref_idx, const sorted_list &query_idx)
{
	task_timer timer("Finding high frequency seeds", 3);
	vector<unsigned> counts(Const::seedp);
	Build_context build_context(ref_idx, query_idx, sid, counts);
	launch_scheduled_thread_pool(build_context, Const::seedp, config.threads_);
	timer.finish();
	log_stream << "Masked positions = " << std::accumulate(counts.begin(), counts.end(), 0) << std::endl;
}
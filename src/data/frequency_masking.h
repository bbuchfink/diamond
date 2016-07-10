/****
Copyright (c) 2014-2015, University of Tuebingen
Author: Benjamin Buchfink
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

#ifndef FREQUENCY_MASKING_H_
#define FREQUENCY_MASKING_H_

#include <algorithm>
#include <vector>
#include "../basic/value.h"

using std::vector;
using std::auto_ptr;

struct Masked_sequence_set : public Sequence_set
{

	Masked_sequence_set(Input_stream &file):
			Sequence_set (file)
	{ }

	Masked_sequence_set()
	{}
	
	void build_masking(unsigned sid, const seedp_range &range, sorted_list &idx)
	{
		task_timer timer ("Counting low complexity seeds", 3);
		vector<unsigned> counts (Const::seedp);
		Count_context count_context (idx, counts);
		launch_scheduled_thread_pool(count_context, Const::seedp, config.threads_);

		timer.finish();
		size_t n = 0;
		for(size_t i=range.begin();i<range.end();++i) {
			n += counts[i];
			const size_t ht_size (std::max(static_cast<size_t>(static_cast<float>(counts[i]) * 1.3), static_cast<size_t>(counts[i] + 1)));
			pos_filters[sid][i] = auto_ptr<filter_table> (new filter_table(ht_size));
		}
		log_stream << "Hit cap = " << config.hit_cap << std::endl;
		log_stream << "Low complexity seeds = " << n << std::endl;

		timer.go("Building position filter");
		Build_context build_context(idx, sid, counts, *this);
		launch_scheduled_thread_pool(build_context, Const::seedp, config.threads_);
		timer.finish();
		log_stream << "Masked positions = " << std::accumulate(counts.begin(), counts.end(), 0) << std::endl;
	}

	bool get_masking(const Letter *pos, unsigned sid) const
	{
		Packed_seed seed;
		shapes[sid].set_seed(seed, pos);
		const filter_table::entry *e;
		if((e = pos_filters[sid][seed_partition(seed)]->operator [](seed_partition_offset(seed))) != 0) {
			const size_t offset (pos - this->data(0));
			return !position_filter(offset, e->value, seed_partition_offset(seed));
		}
		return false;
	}

	virtual ~Masked_sequence_set()
	{ }

private:

	struct Count_context
	{
		Count_context(const sorted_list &idx, vector<unsigned> &counts):
			idx (idx),
			counts (counts)
		{ }
		void operator()(unsigned thread_id, unsigned seedp) const
		{
			unsigned n = 0;
			sorted_list::const_iterator i = idx.get_partition_cbegin(seedp);
			while(!i.at_end()) {
				if(i.n > config.hit_cap)
					++n;
				++i;
			}
			counts[seedp] = n;
		}
		const sorted_list &idx;
		vector<unsigned> &counts;
	};

	struct Build_context
	{
		Build_context(const sorted_list &idx, unsigned sid, vector<unsigned> &counts, Masked_sequence_set &seqs):
			idx (idx),
			sid (sid),
			counts (counts),
			seqs (seqs)
		{ }
		void operator()(unsigned thread_id, unsigned seedp)
		{
			unsigned n = 0;
			sorted_list::iterator i = idx.get_partition_begin(seedp);
			while(!i.at_end()) {
				if(i.n > config.hit_cap)
					n += seqs.mask_seed_pos(i, sid, seedp);
				++i;
			}
			counts[seedp] = n;
		}
		const sorted_list &idx;
		const unsigned sid;
		vector<unsigned> &counts;
		Masked_sequence_set &seqs;
	};

	unsigned mask_seed_pos(sorted_list::iterator &i, unsigned sid, unsigned p)
	{
		const unsigned treshold (filter_treshold((unsigned)i.n));
		unsigned count (0), k (0);
		for(unsigned j=0;j<i.n;++j)
			if(!position_filter(i[j], treshold, i.key())) {
				mask_seed_pos(i[j]);
				++count;
			} else
				*(i.get(k++)) = *(i.get(j));
		i.get(k)->value = 0;
		pos_filters[sid][p]->insert(i.key(), treshold);
		return count;
	}

	void mask_seed_pos(Loc pos)
	{
		Letter *x = this->data(pos);
		*x = set_critical(*x);
	}

private:

	typedef hash_table<uint32_t, uint8_t, value_compare<uint8_t, 0>, murmur_hash> filter_table;
	auto_ptr<filter_table> pos_filters[Const::max_shapes][Const::seedp];

};

#endif /* FREQUENCY_MASKING_H_ */

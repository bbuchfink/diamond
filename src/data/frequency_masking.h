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

#include <vector>

using std::vector;
using std::auto_ptr;

template<typename _val>
struct Masked_sequence_set : public Sequence_set<_val>
{

	Masked_sequence_set(Input_stream &file):
			Sequence_set<_val> (file)
	{ }

	template<typename _loc>
	void build_masking(unsigned sid, const seedp_range &range, typename sorted_list<_loc>::Type &idx)
	{
		task_timer timer ("Counting low complexity seeds", false);
		vector<unsigned> counts (Const::seedp);
		Count_context<_loc> count_context (idx, counts);
		launch_scheduled_thread_pool(count_context, Const::seedp, program_options::threads());

		timer.finish();
		size_t n = 0;
		for(unsigned i=range.begin();i<range.end();++i) {
			n += counts[i];
			const size_t ht_size (std::max(static_cast<size_t>(static_cast<float>(counts[i]) * 1.3), static_cast<size_t>(counts[i] + 1)));
			pos_filters[sid][i] = auto_ptr<filter_table> (new filter_table(ht_size));
		}
		log_stream << "Hit cap = " << program_options::hit_cap << std::endl;
		log_stream << "Low complexity seeds = " << n << std::endl;

		timer.go("Building position filter");
		Build_context<_loc> build_context(idx, sid, counts, *this);
		launch_scheduled_thread_pool(build_context, Const::seedp, program_options::threads());
		timer.finish();
		log_stream << "Masked positions = " << std::accumulate(counts.begin(), counts.end(), 0) << std::endl;
	}

	bool get_masking(const _val *pos, unsigned sid) const
	{
		uint64_t seed;
		shape_config::get().get_shape(sid).set_seed(seed, pos);
		const filter_table::entry *e;
		if((e = pos_filters[sid][seed_partition(seed)]->operator [](seed_partition_offset(seed))) != 0) {
			const size_t offset (pos - this->data(0));
			return !position_filter(offset, e->value, seed_partition_offset(seed));
		} else
			return false;
	}

	virtual ~Masked_sequence_set()
	{ }

private:

	template<typename _loc>
	struct Count_context
	{
		Count_context(const typename sorted_list<_loc>::Type &idx, vector<unsigned> &counts):
			idx (idx),
			counts (counts)
		{ }
		void operator()(unsigned thread_id, unsigned seedp) const
		{
			unsigned n = 0;
			typename sorted_list<_loc>::Type::const_iterator i = idx.get_partition_cbegin(seedp);
			while(!i.at_end()) {
				if(i.n > program_options::hit_cap)
					++n;
				++i;
			}
			counts[seedp] = n;
		}
		const typename sorted_list<_loc>::Type &idx;
		vector<unsigned> &counts;
	};

	template<typename _loc>
	struct Build_context
	{
		Build_context(const typename sorted_list<_loc>::Type &idx, unsigned sid, vector<unsigned> &counts, Masked_sequence_set<_val> &seqs):
			idx (idx),
			sid (sid),
			counts (counts),
			seqs (seqs)
		{ }
		void operator()(unsigned thread_id, unsigned seedp)
		{
			unsigned n = 0;
			typename sorted_list<_loc>::Type::iterator i = idx.get_partition_begin(seedp);
			while(!i.at_end()) {
				if(i.n > program_options::hit_cap)
					n += seqs.mask_seed_pos<_loc>(i, sid, seedp);
				++i;
			}
			counts[seedp] = n;
		}
		const typename sorted_list<_loc>::Type &idx;
		const unsigned sid;
		vector<unsigned> &counts;
		Masked_sequence_set<_val> &seqs;
	};

	template<typename _loc>
	unsigned mask_seed_pos(typename sorted_list<_loc>::Type::iterator &i, unsigned sid, unsigned p)
	{
		const unsigned treshold (filter_treshold(i.n));
		unsigned count (0), k (0);
		for(unsigned j=0;j<i.n;++j)
			if(!position_filter(i[j], treshold, i.key())) {
				mask_seed_pos<_loc>(i[j]);
				++count;
			} else
				*(i.get(k++)) = *(i.get(j));
		i.get(k)->value = 0;
		pos_filters[sid][p]->insert(i.key(), treshold);
		return count;
	}

	template<typename _loc>
	void mask_seed_pos(_loc pos)
	{
		_val *x = this->data(pos);
		*x = set_critical(*x);
	}

private:

	typedef hash_table<uint32_t, uint8_t, value_compare<uint8_t, 0>, murmur_hash> filter_table;
	auto_ptr<filter_table> pos_filters[Const::max_shapes][Const::seedp];

};

#endif /* FREQUENCY_MASKING_H_ */

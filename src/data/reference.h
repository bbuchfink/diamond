/****
Copyright (c) 2014, University of Tuebingen
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****
Author: Benjamin Buchfink
****/

#ifndef REFERENCE_H_
#define REFERENCE_H_

#include <memory>
#include <string>
#include <numeric>
#include "../util/binary_file.h"
#include "sorted_list.h"
#include "../basic/statistics.h"
#include "../data/seed_histogram.h"
#include "../util/hash_table.h"
#include "../util/hash_function.h"
#include "../basic/packed_loc.h"
#include "sequence_set.h"
#include "boost/ptr_container/ptr_vector.hpp"

using std::auto_ptr;
using boost::ptr_vector;

struct Reference_header
{
	Reference_header():
		unique_id (0x24af8a415ee186dllu),
		build (Const::build_version),
		long_addressing (false),
		sequences (0),
		letters (0)
	{ }
	uint64_t unique_id;
	uint32_t build;
	bool long_addressing;
	unsigned n_blocks;
	size_t sequences, letters;
	double block_size;
#ifdef EXTRA
	Sequence_type sequence_type;
#endif
} ref_header;

struct Database_format_exception : public exception
{
	virtual const char* what() const throw()
	{ return "Database file is not a DIAMOND database."; }
};

template<typename _val>
struct Database_file : public Input_stream
{
	Database_file():
		Input_stream (program_options::database)
	{
		if(this->read(&ref_header, 1) != 1)
			throw Database_format_exception ();
		if(ref_header.unique_id != Reference_header ().unique_id)
			throw Database_format_exception ();
		if(ref_header.build > Const::build_version || ref_header.build < Const::build_compatibility)
			throw invalid_database_version_exception();
#ifdef EXTRA
		if(sequence_type(_val()) != ref_header.sequence_type)
			throw std::runtime_error("Database has incorrect sequence type for current alignment mode.");
#endif
	}
	void rewind()
	{ this->seek(sizeof(Reference_header)); }
};

template<typename _val>
struct ref_seqs
{
	static const Sequence_set<_val>& get()
	{ return *data_; }
	static Sequence_set<_val> *data_;
};

template<typename _val> Sequence_set<_val>* ref_seqs<_val>::data_ = 0;

struct ref_ids
{
	static const String_set<char,0>& get()
	{ return *data_; }
	static String_set<char,0> *data_;
};

String_set<char,0>* ref_ids::data_ = 0;

seed_histogram ref_hst;
unsigned current_ref_block;

size_t max_id_len(const String_set<char,0> &ids)
{
	size_t max (0);
	for(size_t i=0;i<ids.get_length(); ++i)
		max = std::max(max, find_first_of(ids[i].c_str(), Const::id_delimiters));
	return max;
}

struct Masking
{

	template<typename _val, typename _loc>
	void build(unsigned sid, const seedp_range &range, typename sorted_list<_loc>::Type &idx)
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
		Build_context<_val,_loc> build_context(idx, sid, counts, *this);
		launch_scheduled_thread_pool(build_context, Const::seedp, program_options::threads());
		timer.finish();
		log_stream << "Masked positions = " << std::accumulate(counts.begin(), counts.end(), 0) << std::endl;
	}

	template<typename _val>
	bool get(const _val *pos, unsigned sid) const
	{
		uint64_t seed;
		shape_config::get().get_shape(sid).set_seed(seed, pos);
		const filter_table::entry *e;
		if((e = pos_filters[sid][seed_partition(seed)]->operator [](seed_partition_offset(seed))) != 0) {
			const size_t offset (pos - ref_seqs<_val>::data_->data(0));
			return !position_filter(offset, e->value, seed_partition_offset(seed));
		} else
			return false;
	}

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

	template<typename _val, typename _loc>
	struct Build_context
	{
		Build_context(const typename sorted_list<_loc>::Type &idx, unsigned sid, vector<unsigned> &counts, Masking &masking):
			idx (idx),
			sid (sid),
			counts (counts),
			masking (masking)
		{ }
		void operator()(unsigned thread_id, unsigned seedp)
		{
			unsigned n = 0;
			typename sorted_list<_loc>::Type::iterator i = idx.get_partition_begin(seedp);
			while(!i.at_end()) {
				if(i.n > program_options::hit_cap)
					n += masking.mask_seed_pos<_val,_loc>(i, sid, seedp);
				++i;
			}
			counts[seedp] = n;
		}
		const typename sorted_list<_loc>::Type &idx;
		const unsigned sid;
		vector<unsigned> &counts;
		Masking &masking;
	};

	template<typename _val, typename _loc>
	unsigned mask_seed_pos(typename sorted_list<_loc>::Type::iterator &i, unsigned sid, unsigned p)
	{
		const unsigned treshold (filter_treshold(i.n));
		unsigned count (0), k (0);
		for(unsigned j=0;j<i.n;++j)
			if(!position_filter(i[j], treshold, i.key())) {
				mask_seed_pos<_val,_loc>(i[j]);
				++count;
			} else
				*(i.get(k++)) = *(i.get(j));
		i.get(k)->value = 0;
		pos_filters[sid][p]->insert(i.key(), treshold);
		return count;
	}

	template<typename _val, typename _loc>
	static void mask_seed_pos(_loc pos)
	{
		_val *x (ref_seqs<_val>::data_->data(pos));
		*x = set_critical(*x);
	}

private:

	typedef hash_table<uint32_t, uint8_t, value_compare<uint8_t, 0>, murmur_hash> filter_table;
	auto_ptr<filter_table> pos_filters[Const::max_shapes][Const::seedp];

} ref_masking;

struct Ref_map
{
	Ref_map():
		next_ (0)
	{ }
	void init(unsigned ref_count)
	{
		const unsigned block = current_ref_block;
		if(data_.size() < block+1) {
			data_.resize(block+1);
			data_[block].insert(data_[block].end(), ref_count, std::numeric_limits<unsigned>::max());
		}
	}
	template<typename _val>
	uint32_t get(unsigned block, unsigned i)
	{
		uint32_t n = data_[block][i];
		if(n != std::numeric_limits<unsigned>::max())
			return n;
		else {
			tthread::lock_guard<tthread::mutex> lock (mtx_);
			n = data_[block][i];
			if(n != std::numeric_limits<uint32_t>::max())
				return n;
			n = next_++;
			data_[block][i] = n;
			len_.push_back(ref_seqs<_val>::get().length(i));
			name_.push_back(get_str(ref_ids::get()[i].c_str(), Const::id_delimiters));
			return n;
		}
	}
	/*template<typename _val>
	void finish()
	{
		vector<pair<unsigned,unsigned> > v;
		for(unsigned i=0;i<data_.size();++i)
			if(data_[i] != std::numeric_limits<unsigned>::max())
				v.push_back(pair<unsigned,unsigned> (data_[i], i));
		std::sort(v.begin(), v.end());
		for(vector<pair<unsigned,unsigned> >::const_iterator i = v.begin(); i!=v.end(); ++i) {
			const char* s = ref_ids::get()[i->second].c_str();
			buf_ << (uint32_t)(strlen(s)+1);
			buf_.write_c_str(s);
			buf_ << (uint32_t)ref_seqs<_val>::get()[i->second].length();
		}
	}*/
private:
	tthread::mutex mtx_;
	vector<vector<uint32_t> > data_;
	vector<uint32_t> len_;
	ptr_vector<string> name_;
	uint32_t next_;
	friend struct DAA_output;
} ref_map;

#endif /* REFERENCE_H_ */

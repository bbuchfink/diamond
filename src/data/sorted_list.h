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

#ifndef SORTED_LIST_H_
#define SORTED_LIST_H_

#include "../util/util.h"
#include "seed_histogram.h"
#include "../basic/packed_loc.h"

template<typename _pos>
struct sorted_list
{

	typedef sorted_list<typename packed_sequence_location<_pos>::type> Type;

	struct entry
	{
		entry():
			key (),
			value ()
		{ }
		entry(unsigned key, _pos value):
			key (key),
			value (value)
		{ }
		bool operator<(const entry &rhs) const
		{ return key < rhs.key; }
		unsigned	key;
		_pos		value;
	} __attribute__((packed));

	static char* alloc_buffer(const seed_histogram &hst)
	{ return new char[sizeof(entry) * hst.max_chunk_size()]; }

	template<typename _val>
	sorted_list(char *buffer, const Sequence_set<_val> &seqs, const shape &sh, const shape_histogram &hst, const seedp_range &range):
		limits_ (hst, range),
		data_ (reinterpret_cast<entry*>(buffer))
	{
		task_timer timer ("Building seed list", false);
		Build_context<_val> build_context (seqs, sh, range, build_iterators(hst));
		launch_scheduled_thread_pool(build_context, Const::seqp, program_options::threads());
		timer.go("Sorting seed list");
		Sort_context sort_context (*this);
		launch_scheduled_thread_pool(sort_context, Const::seedp, program_options::threads());
	}

	template<typename _t>
	struct Iterator_base
	{
		Iterator_base(_t *i, _t *end):
			i (i),
			end (end),
			n (count())
		{ }
		size_t count() const
		{
			_t *k (i);
			size_t n (0);
			while(k < end && (k++)->key == i->key)
				++n;
			return n;
		}
		void operator++()
		{ i += n; n = count(); }
		_pos operator[](unsigned k) const
		{ return (i+k)->value; }
		_t* get(unsigned k)
		{ return i+k; }
		bool at_end() const
		{ return i >= end; }
		unsigned key() const
		{ return i->key; }
		_t *i, *end;
		size_t n;
	};

	typedef Iterator_base<entry> iterator;
	typedef Iterator_base<const entry> const_iterator;

	const_iterator get_partition_cbegin(unsigned p) const
	{ return const_iterator (cptr_begin(p), cptr_end(p)); }

	iterator get_partition_begin(unsigned p) const
	{ return iterator (ptr_begin(p), ptr_end(p)); }

private:

	typedef Static_matrix<entry*,Const::seqp,Const::seedp> Ptr_set;

	struct buffered_iterator
	{
		static const unsigned BUFFER_SIZE = 16;
		buffered_iterator(entry **ptr)
		{
			memset(n, 0, sizeof(n));
			memcpy(this->ptr, ptr, sizeof(this->ptr));
		}
		void push(seed key, _pos value, const seedp_range &range)
		{
			const unsigned p (seed_partition(key));
			if(range.contains(p)) {
				assert(n[p] < BUFFER_SIZE);
				buf[p][n[p]++] = entry (seed_partition_offset(key), value);
				if(n[p] == BUFFER_SIZE)
					flush(p);
			}
		}
		void flush(unsigned p)
		{
			memcpy(ptr[p], buf[p], n[p] * sizeof(entry));
			ptr[p] += n[p];
			n[p] = 0;
		}
		void flush()
		{
			for(unsigned p=0;p<Const::seedp;++p)
				if(n[p] > 0)
					flush(p);
		}
		entry* ptr[Const::seedp];
		entry  	 buf[Const::seedp][BUFFER_SIZE];
		uint8_t  n[Const::seedp];
	};

	entry* ptr_begin(unsigned i) const
	{ return &data_[limits_[i]]; }

	entry* ptr_end(unsigned i) const
	{ return &data_[limits_[i+1]]; }

	const entry* cptr_begin(unsigned i) const
	{ return &data_[limits_[i]]; }

	const entry* cptr_end(unsigned i) const
	{ return &data_[limits_[i+1]]; }

	template<typename _val>
	struct Build_context
	{
		Build_context(const Sequence_set<_val> &seqs, const shape &sh, const seedp_range &range, Ptr_set *iterators):
			seqs (seqs),
			sh (sh),
			range (range),
			iterators (iterators),
			seq_partition (seqs.partition())
		{ }
		void operator()(unsigned thread_id, unsigned seqp) const
		{
			build_seqp<_val>(seqs,
					seq_partition[seqp],
					seq_partition[seqp+1],
					(*iterators)[seqp],
					sh,
					range);
		}
		const Sequence_set<_val> &seqs;
		const shape &sh;
		const seedp_range &range;
		const auto_ptr<Ptr_set> iterators;
		const vector<size_t> seq_partition;
	};

	template<typename _val>
	static void build_seqp(const Sequence_set<_val> &seqs, unsigned begin, unsigned end, entry **ptr, const shape &sh, const seedp_range &range)
	{
		uint64_t key;
		auto_ptr<buffered_iterator> it (new buffered_iterator(ptr));
		for(size_t i=begin;i<end;++i) {
			const sequence<const _val> seq = seqs[i];
			if(seq.length()<sh.length_) continue;
			for(unsigned j=0;j<seq.length()-sh.length_+1; ++j) {
				if(sh.set_seed(key, &seq[j]))
					it->push(key, seqs.position(i, j), range);
			}
		}
		it->flush();
	}

	Ptr_set* build_iterators(const shape_histogram &hst) const
	{
		Ptr_set *iterators = new Ptr_set;
		for(unsigned i=0;i<Const::seedp;++i)
			(*iterators)[0][i] = ptr_begin(i);

		for(unsigned i=1;i<Const::seqp;++i) {
			for(unsigned j=0;j<Const::seedp;++j)
				(*iterators)[i][j] = (*iterators)[i-1][j] + hst[i-1][j];
		}
		return iterators;
	}

	struct Sort_context
	{
		Sort_context(sorted_list &sl):
			sl (sl)
		{ }
		void operator()(unsigned thread_id ,unsigned seedp) const
		{ std::sort(sl.ptr_begin(seedp), sl.ptr_end(seedp)); }
		sorted_list &sl;
	};

	struct Limits : vector<size_t>
	{
		Limits(const shape_histogram &hst, const seedp_range &range)
		{
			task_timer timer ("Computing limits", false);
			this->push_back(0);
			for(unsigned i=0;i<Const::seedp;++i)
				this->push_back(this->operator[](i) + (range.contains(i) ? partition_size(hst, i) : 0));
		}
	};

	const Limits limits_;
	entry *data_;

};

#endif /* SORTED_LIST_H_ */

/****
Copyright (c) 2014-2016, University of Tuebingen, Benjamin Buchfink
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

#ifndef SORTED_LIST_H_
#define SORTED_LIST_H_

#include "../util/util.h"
#include "seed_histogram.h"
#include "../basic/packed_loc.h"
#include "../util/system.h"

#pragma pack(1)

struct sorted_list
{

	typedef Packed_loc _pos;

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
	} PACKED_ATTRIBUTE ;

	static char* alloc_buffer(const Partitioned_histogram &hst)
	{ return new char[sizeof(entry) * hst.max_chunk_size()]; }

	sorted_list()
	{}

	sorted_list(char *buffer, const Sequence_set &seqs, const shape &sh, const shape_histogram &hst, const seedp_range &range, const vector<size_t> seq_partition):
		limits_ (hst, range),
		data_ (reinterpret_cast<entry*>(buffer))
	{
		task_timer timer ("Building seed list", 3);
		Build_context build_context (seqs, sh, range, build_iterators(hst), seq_partition);
		launch_scheduled_thread_pool(build_context, unsigned(seq_partition.size() - 1), config.threads_);
		timer.go("Sorting seed list");
		Sort_context sort_context (*this);
		launch_scheduled_thread_pool(sort_context, Const::seedp, config.threads_);
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
		Loc operator[](unsigned k) const
		{ return (Loc)((i+k)->value); }
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

	struct Random_access_iterator
	{
		Random_access_iterator(const entry *i, const entry *end) :
			i(i),
			end(end),
			key_(i ? i->key : 0)
		{ }
		void operator++()
		{
			++i;
		}
		Loc operator*() const
		{
			return (Loc)i->value;
		}
		bool good() const
		{
			return i < end && i->key == key_;
		}
		unsigned key() const
		{
			return key_;
		}
		const entry *i, *end;
		unsigned key_;
	};

	size_t iterator_offset(const const_iterator &i, unsigned p) const
	{
		return i.i - cptr_begin(p);
	}

	Random_access_iterator random_access(unsigned p, size_t offset) const
	{
		return Random_access_iterator(cptr_begin(p) + offset, cptr_end(p));
	}

private:

	typedef vector<Array<entry*, Const::seedp> > Ptr_set;

	struct buffered_iterator
	{
		static const unsigned BUFFER_SIZE = 16;
		buffered_iterator(entry* const* ptr)
		{
			memset(n, 0, sizeof(n));
			memcpy(this->ptr, ptr, sizeof(this->ptr));
		}
		void push(Packed_seed key, Loc value, const seedp_range &range)
		{
			const unsigned p(seed_partition(key));
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

	struct Build_context
	{
		Build_context(const Sequence_set &seqs, const shape &sh, const seedp_range &range, const Ptr_set &iterators, const vector<size_t> &seq_partition):
			seqs (seqs),
			sh (sh),
			range (range),
			iterators (iterators),
			seq_partition (seq_partition)
		{ }
		void operator()(unsigned thread_id, unsigned seqp) const
		{
			build_seqp(seqs,
				seq_partition[seqp],
				seq_partition[seqp + 1],
				iterators[seqp].begin(),
				sh,
				range);
		}
		const Sequence_set &seqs;
		const shape &sh;
		const seedp_range &range;
		const Ptr_set iterators;
		const vector<size_t> &seq_partition;
	};

	static void build_seqp(const Sequence_set &seqs, size_t begin, size_t end, entry* const* ptr, const shape &sh, const seedp_range &range)
	{
		uint64_t key;
		auto_ptr<buffered_iterator> it (new buffered_iterator(ptr));
		for(size_t i=begin;i<end;++i) {
			const sequence seq = seqs[i];
			if(seq.length()<sh.length_) continue;
			for (unsigned j = 0; j < seq.length() - sh.length_ + 1; ++j) {
				if (sh.set_seed(key, &seq[j]))
					it->push(key, seqs.position(i, j), range);
			}
		}
		it->flush();
	}

	Ptr_set build_iterators(const shape_histogram &hst) const
	{
		Ptr_set iterators (hst.size());
		for (unsigned i = 0; i < Const::seedp; ++i)
			iterators[0][i] = ptr_begin(i);

		for (unsigned i = 1; i < hst.size(); ++i)
			for (unsigned j = 0; j < Const::seedp; ++j)
				iterators[i][j] = iterators[i - 1][j] + hst[i - 1][j];
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
		Limits()
		{}
		Limits(const shape_histogram &hst, const seedp_range &range)
		{
			task_timer timer ("Computing limits", 3);
			this->push_back(0);
			for(unsigned i=0;i<Const::seedp;++i) {
#ifdef EXTRA
				log_stream << i << ' ' << partition_size(hst, i) << endl;
#endif
				this->push_back(this->operator[](i) + (range.contains(i) ? partition_size(hst, i) : 0));
			}
		}
	};

	Limits limits_;
	entry *data_;

};

template<typename _it>
struct Merge_iterator
{
	Merge_iterator(const _it &i, const _it& j):
		i(i),
		j(j)
	{}
	bool next()
	{		
		while (!i.at_end() && !j.at_end()) {
			if (i.key() < j.key()) {
				++i;
			}
			else if (j.key() < i.key()) {
				++j;
			}
			else
				return true;
		}
		return false;
	}
	void operator++()
	{
		++i;
		++j;
	}
	_it i, j;
};

#pragma pack()

#endif /* SORTED_LIST_H_ */

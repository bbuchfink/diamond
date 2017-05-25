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

	static char* alloc_buffer(const Partitioned_histogram &hst);
	sorted_list();
	sorted_list(char *buffer, const Sequence_set &seqs, const shape &sh, const shape_histogram &hst, const seedp_range &range, const vector<size_t> seq_partition, const Seed_set *filter = 0);

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

	const_iterator get_partition_cbegin(unsigned p) const;
	iterator get_partition_begin(unsigned p) const;

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

	Random_access_iterator random_access(unsigned p, size_t offset) const;

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

	entry* ptr_begin(unsigned i) const;
	entry* ptr_end(unsigned i) const;
	const entry* cptr_begin(unsigned i) const;
	const entry* cptr_end(unsigned i) const;

	struct Build_context
	{
		Build_context(const Sequence_set &seqs, const shape &sh, const seedp_range &range, const Ptr_set &iterators, const vector<size_t> &seq_partition, const Seed_set *filter):
			seqs (seqs),
			sh (sh),
			range (range),
			iterators (iterators),
			seq_partition (seq_partition),
			filter(filter)
		{ }
		void operator()(unsigned thread_id, unsigned seqp) const
		{
			build_seqp(seqs,
				seq_partition[seqp],
				seq_partition[seqp + 1],
				iterators[seqp].begin(),
				sh,
				range,
				filter);
		}
		const Sequence_set &seqs;
		const shape &sh;
		const seedp_range &range;
		const Ptr_set iterators;
		const vector<size_t> &seq_partition;
		const Seed_set *filter;
	};

	static void build_seqp(const Sequence_set &seqs, size_t begin, size_t end, entry* const* ptr, const shape &sh, const seedp_range &range, const Seed_set *filter);
	Ptr_set build_iterators(const shape_histogram &hst) const;

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

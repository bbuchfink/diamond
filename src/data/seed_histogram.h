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

#ifndef SEED_HISTOGRAM_H_
#define SEED_HISTOGRAM_H_

#include <limits>
#include "../basic/seed.h"
#include "sequence_set.h"
#include "../basic/shape_config.h"
#include "../util/thread.h"
#include "../basic/seed_iterator.h"
#include "seed_set.h"

using std::vector;

typedef vector<Array<unsigned,Const::seedp> > shape_histogram;

struct seedp_range
{
	seedp_range():
		begin_ (0),
		end_ (0)
	{ }
	seedp_range(unsigned begin, unsigned end):
		begin_ (begin),
		end_ (end)
	{ }
	bool contains(unsigned i) const
	{ return i >= begin_ && i < end_; }
	unsigned begin() const
	{ return begin_; }
	unsigned end() const
	{ return end_; }
	bool lower(unsigned i) const
	{ return i < begin_; }
	bool lower_or_equal(unsigned i) const
	{ return i < end_; }
	unsigned size() const
	{
		return end_ - begin_;
	}
	static seedp_range all()
	{
		return seedp_range(0, Const::seedp);
	}
private:
	unsigned begin_, end_;
};

extern seedp_range current_range;

inline size_t partition_size(const shape_histogram &hst, unsigned p)
{
	size_t s = 0;
	for(unsigned i=0;i<hst.size();++i)
		s += hst[i][p];
	return s;
}

inline size_t hst_size(const shape_histogram &hst, const seedp_range &range)
{
	size_t s = 0;
	for(unsigned i=range.begin();i<range.end();++i)
		s += partition_size(hst, i);
	return s;
}

struct Partitioned_histogram
{

	Partitioned_histogram();
	
	template<typename _filter>
	Partitioned_histogram(const Sequence_set &seqs, bool serial, const _filter *filter) :
		data_(shapes.count()),
		p_(seqs.partition(config.threads_))
	{
		for (unsigned s = 0; s < shapes.count(); ++s) {
			data_[s].resize(p_.size() - 1);
			memset(data_[s].data(), 0, (p_.size() - 1)*sizeof(unsigned)*Const::seedp);
		}
		Ptr_vector<Callback> cb;
		for (size_t i = 0; i < p_.size() - 1; ++i)
			cb.push_back(new Callback(i, data_));
		if (serial)
			for (unsigned s = 0; s < shapes.count(); ++s)
				seqs.enum_seeds(cb, p_, s, s + 1, filter);
		else
			seqs.enum_seeds(cb, p_, 0, shapes.count(), filter);
	}

	const shape_histogram& get(unsigned sid) const
	{ return data_[sid]; }

	size_t max_chunk_size() const;

	const vector<size_t>& partition() const
	{
		return p_;
	}

private:

	struct Callback
	{
		Callback(size_t seqp, vector<shape_histogram> &data)
		{
			for (unsigned s = 0; s < shapes.count(); ++s)
				ptr.push_back(data[s][seqp].begin());
		}
		bool operator()(uint64_t seed, uint64_t pos, size_t shape)
		{
			++ptr[shape][seed_partition(seed)];
			return true;
		}
		void finish() const
		{
		}
		vector<unsigned*> ptr;
	};

	vector<shape_histogram> data_;
	vector<size_t> p_;

};

#endif /* SEED_HISTOGRAM_H_ */
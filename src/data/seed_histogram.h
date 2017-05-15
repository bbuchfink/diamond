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

	Partitioned_histogram()
	{ }

	Partitioned_histogram(const Sequence_set &seqs, unsigned longest):
		data_(shapes.count()),
		p_(seqs.partition(config.threads_*4))
	{
		for (unsigned s = 0; s < shapes.count(); ++s) {
			data_[s].resize(p_.size() - 1);
			memset(data_[s].data(), 0, (p_.size() - 1)*sizeof(unsigned)*Const::seedp);
		}
		Build_context context (seqs, *this, longest);
		launch_scheduled_thread_pool(context, (unsigned)(p_.size() - 1), (unsigned)(p_.size() - 1));
	}

	const shape_histogram& get(unsigned sid) const
	{ return data_[sid]; }

	size_t max_chunk_size() const
	{
		size_t max = 0;
		::partition<unsigned> p(Const::seedp, config.lowmem);
		for (unsigned shape = 0; shape < shapes.count(); ++shape)
			for (unsigned chunk = 0; chunk < p.parts; ++chunk)
				max = std::max(max, hst_size(data_[shape], seedp_range(p.getMin(chunk), p.getMax(chunk))));
		return max;
	}

	const vector<size_t>& partition() const
	{
		return p_;
	}

private:

	struct Build_context
	{
		Build_context(const Sequence_set &seqs, Partitioned_histogram &hst, unsigned longest):
			seqs (seqs),
			hst (hst),
			longest(longest)
		{ }
		void operator()(unsigned thread_id, unsigned seqp) const
		{
			vector<char> buf(longest);
			hst.build_seq_partition(seqs, seqp, hst.p_[seqp], hst.p_[seqp+1], buf);
		}
		const Sequence_set &seqs;
		Partitioned_histogram &hst;
		unsigned longest;
	};

	void build_seq_partition(const Sequence_set &seqs,
		const unsigned seqp,
		const size_t begin,
		const size_t end,
		vector<char> &buf)
	{
		for (size_t i = begin; i < end; ++i) {

			assert(i < seqs.get_length());
			if (seqs[i].length() == 0)
				continue;
			Reduction::reduce_seq(seqs[i], buf);

			for (unsigned s = 0; s < shapes.count(); ++s) {
				Seed_iterator it(buf, shapes[s]);
				unsigned *ptr = data_[s][seqp].begin();
				uint64_t seed;
				while (it.good())
					if (it.get(seed, shapes[s]))
						++ptr[seed_partition(seed)];
			}

		}
	}

	vector<shape_histogram> data_;
	vector<size_t> p_;

};

#endif /* SEED_HISTOGRAM_H_ */
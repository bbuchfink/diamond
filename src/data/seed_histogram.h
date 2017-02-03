/****
Copyright (c) 2014-16, University of Tuebingen, Benjamin Buchfink
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
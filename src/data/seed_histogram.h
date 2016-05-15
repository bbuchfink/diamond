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

#ifndef SEED_HISTOGRAM_H_
#define SEED_HISTOGRAM_H_

#include <limits>
#include "../basic/seed.h"
#include "sequence_set.h"
#include "../basic/shape_config.h"
#include "../util/thread.h"

using std::vector;

inline void encode_zero_rle(const int32_t *data, size_t len, Output_stream &out)
{
	const int32_t *p = data, *end = data + len;
	int32_t n = 0;
	while(p < end) {
		while(p < end && *p == 0) {
			--n;
			++p;
		}
		if(n < 0) {
			out.write(&n, 1);
			n = 0;
		}
		if(p < end) {
			out.write(p, 1);
			++p;
		}
	}
}

inline void decode_zero_rle(int32_t *data, size_t len, Input_stream &in)
{
	size_t n = 0;
	while(n < len) {
		int32_t x;
		in.read(&x, 1);
		if(x >= 0) {
			*(data++) = x;
			++n;
		} else {
			for(int32_t i=0;i>x;--i) {
				*(data++) = 0;
				++n;
			}
		}
	}
}

typedef int32_t shape_histogram[Const::seqp][Const::seedp];

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
private:
	unsigned begin_, end_;
};

extern seedp_range current_range;

inline size_t partition_size(const shape_histogram &hst, unsigned p)
{
	size_t s (0);
	for(unsigned i=0;i<Const::seqp;++i)
		s += hst[i][p];
	return s;
}

inline size_t hst_size(const shape_histogram &hst, const seedp_range &range)
{
	size_t s (0);
	for(unsigned i=range.begin();i<range.end();++i)
		s += partition_size(hst, i);
	return s;
}

struct seed_histogram
{

	seed_histogram()
	{ }

	seed_histogram(const Sequence_set &seqs)
	{
		memset(data_, 0, sizeof(data_));
		Build_context context (seqs, *this);
		launch_scheduled_thread_pool(context, Const::seqp, config.threads_);
	}

	const shape_histogram& get(unsigned index_mode, unsigned sid) const
	{ return data_[index_mode][sid]; }

	size_t max_chunk_size() const
	{
		size_t max (0);
		::partition<unsigned> p (Const::seedp, config.lowmem);
		for(unsigned shape=0;shape < shapes.count();++shape)
			for(unsigned chunk=0;chunk < p.parts; ++chunk)
				max = std::max(max, hst_size(data_[config.index_mode][shape], seedp_range(p.getMin(chunk), p.getMax(chunk))));
		return max;
	}

	void save(Output_stream &out) const
	{
		encode_zero_rle(reinterpret_cast<const int32_t*>(data_), sizeof(data_)/sizeof(int32_t), out);
	}

	void load(Input_stream &in)
	{
		decode_zero_rle(reinterpret_cast<int32_t*>(data_), sizeof(data_)/sizeof(int32_t), in);
	}

private:

	struct Build_context
	{
		Build_context(const Sequence_set &seqs, seed_histogram &hst):
			seqs (seqs),
			cfgs (shape_configs()),
			seq_partition (seqs.partition()),
			hst (hst)
		{ }
		void operator()(unsigned thread_id, unsigned seqp) const
		{ hst.build_seq_partition(seqs, seqp, seq_partition[seqp], seq_partition[seqp+1], cfgs); }
		const Sequence_set &seqs;
		const vector<shape_config> cfgs;
		const vector<size_t> seq_partition;
		seed_histogram &hst;
	};

	void build_seq_partition(const Sequence_set &seqs,
			const unsigned seqp,
			const size_t begin,
			const size_t end,
			const vector<shape_config> &cfgs)
	{
		assert(seqp < Const::seqp);
		uint64_t key;
		for(size_t i=begin;i<end;++i) {

			assert(i < seqs.get_length());
			const sequence seq = seqs[i];
			if(seq.length() < Const::min_shape_len) continue;
			for(unsigned j=0;j<seq.length()+1-Const::min_shape_len; ++j)
				for(vector<shape_config>::const_iterator cfg = cfgs.begin(); cfg != cfgs.end(); ++cfg) {
					assert(cfg->mode() < Const::index_modes);
					assert(cfg->count() <= Const::max_shapes);
					for(unsigned k=0;k<cfg->count(); ++k)
						if(j+cfg->get_shape(k).length_ < seq.length()+1 && cfg->get_shape(k).set_seed(key, &seq[j]))
							++data_[cfg->mode()][k][seqp][seed_partition(key)];
				}

		}
	}

	static vector<shape_config> shape_configs()
	{
		vector<shape_config> v;
		if(config.command == Config::makedb) {
			for(unsigned i=0;i<Const::index_modes;++i)
				v.push_back(shape_config (i, Const::max_shapes));
		} else
			v.push_back(shape_config (config.index_mode, Const::max_shapes));
		return v;
	}

	shape_histogram data_[Const::index_modes][Const::max_shapes];

};

#endif /* SEED_HISTOGRAM_H_ */

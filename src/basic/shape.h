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

#ifndef SHAPE_H_
#define SHAPE_H_

#include <string.h>
#include "const.h"
#include "value.h"
#include "seed.h"
#include "score_matrix.h"
#include "reduction.h"
#include "../util/util.h"
#include "config.h"

/*struct All_partitions {};
struct Filter_partition {};

template<typename _f>
bool include_partition(unsigned p)
{
	return true;
}

template<>
bool include_partition<Filter_partition>(unsigned p)
{
	return current_range.contains(p);
}*/


struct Letter_trail
{
	Letter_trail()
	{
		bucket[0] = 0;
		for (int i = 1; i < 20; ++i)
			bucket[i] = -1;
	}
	Letter_trail(const Reduction &reduction)
	{
		for (int i = 0; i < 20; ++i)
			bucket[i] = reduction(i);
	}
	int operator()(char l) const
	{
		return bucket[(long)l];
	}
	int next() const
	{
		for (int i = 0; i < 20; ++i)
			if (bucket[i] == -1)
				return i;
		return -1;
	}
	int buckets() const
	{
		int m = 0;
		for (int i = 0; i < 20; ++i)
			m = std::max(m, bucket[i]);
		return m + 1;
	}
	double background_p() const;
	double foreground_p(double id) const;
	friend std::ostream& operator<<(std::ostream &s, const Letter_trail &t);
	int bucket[20];
};

#define OPT_W 7
typedef Letter_trail Trail[OPT_W];

struct shape
{

	shape():
		length_ (0),
		weight_ (0),
		d_ (0),
		mask_ (0),
		rev_mask_ (0),
		id_ (0)
	{ memset(positions_, 0, sizeof(uint32_t)*Const::max_seed_weight); }

	shape(const char *code, unsigned id):
		weight_ (0),
		mask_ (0),
		rev_mask_ (0),
		id_ (id)
	{
		assert(id < Const::max_shapes);
		assert(strlen(code) <= 32);
		memset(positions_, 0, sizeof(uint32_t)*Const::max_seed_weight);
		unsigned i (0);
		for(;i<strlen(code);++i) {
			rev_mask_ <<= 1;
			if(code[i] == '1') {
				assert(weight_ < Const::max_seed_weight);
				positions_[weight_] = i;
				++weight_;
				mask_ |= 1 << i;
				rev_mask_ |= 1;
			}
		}
		length_ = i;
		d_ = positions_[weight_/2-1];
	}

	inline bool set_seed(Packed_seed &s, const Letter *seq) const
	{
		s = 0;
#ifdef FREQUENCY_MASKING
		double f = 0;
#endif
		for(unsigned i=0;i<weight_;++i) {
			Letter l = seq[positions_[i]];
			if (l == value_traits.mask_char || l == '\xff')
				return false;
			unsigned r = Reduction::reduction(l);
#ifdef FREQUENCY_MASKING
			f += background_freq[r];
#endif
			s *= Reduction::reduction.size();
			s += uint64_t(r);
		}
#ifdef FREQUENCY_MASKING
		if(use_seed_freq() && f > config.max_seed_freq) return false;
#endif
		return true;
	}

	inline bool set_seed_reduced(Packed_seed &s, const Letter *seq) const
	{
		s = 0;
		for (unsigned i = 0; i < weight_; ++i) {
			Letter l = seq[positions_[i]];
			if (l == value_traits.mask_char)
				return false;
			s *= Reduction::reduction.size();
			s += uint64_t(l);
		}
		return true;
	}

	inline bool set_seed(Seed &s, const Letter *seq) const
	{
		for (unsigned i = 0; i < weight_; ++i) {
			Letter l = seq[positions_[i]];
			if (l >= 20)
				return false;
			s[i] = l;
		}
		return true;
	}
	
	friend std::ostream& operator<<(std::ostream&s, const shape &sh)
	{
		for (unsigned i = 0; i < sh.length_; ++i)
			s << ((sh.mask_ & (1 << i)) ? '1' : '0');
		return s;
	}

	int score(const Letter *x, const Letter *y) const
	{
		int score = 0;
		for (unsigned i = 0; i < weight_; ++i)
			score += score_matrix(x[positions_[i]], y[positions_[i]]);
		return score;
	}

	bool hit(const Letter *x, const Letter *y, const Trail &trail) const
	{
		for (unsigned i = 0; i < weight_; ++i)
#if OPT_W==1
			if (trail[0](x[positions_[i]]) != trail[0](y[positions_[i]]))
				return false;
#else
			if (trail[i](x[positions_[i]]) != trail[i](y[positions_[i]]))
				return false;
#endif
		return true;
	}

	uint32_t length_, weight_, positions_[Const::max_seed_weight], d_, mask_, rev_mask_, id_;

};

#endif /* SHAPE_H_ */
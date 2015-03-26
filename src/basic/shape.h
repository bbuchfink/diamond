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

#ifndef SHAPE_H_
#define SHAPE_H_

#include "const.h"
#include "value.h"
#include "seed.h"
#include "../util/hash_function.h"
#include "score_matrix.h"
#include "reduction.h"

const double background_freq[] = {-1.188861,
		-4.343446,
		-2.648093,
		-3.806941,
		-3.742636,
		-3.221182,
		-3.498273,
		-1.498637,
		-4.339607,
		-3.027002,
		-1.557546 };

template<typename _val>
bool use_seed_freq()
{ return false; }

template<>
bool use_seed_freq<Amino_acid>()
{ return true; }

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

	template<typename _val>
	inline bool set_seed(uint64_t &s, const _val *seq) const
	{
		s = 0;
		double f = 0;
		for(unsigned i=0;i<weight_;++i) {
			_val l = seq[positions_[i]];
			if(l == Value_traits<_val>::MASK_CHAR || l == String_set<_val>::PADDING_CHAR)
				return false;
			l = mask_critical(l);
			unsigned r = Reduction<_val>::reduction(l);
			f += background_freq[r];
			s *= Reduction<_val>::reduction.size();
			s += uint64_t(r);
		}
		if(use_seed_freq<_val>() && f > program_options::max_seed_freq) return false;
		return true;
	}

	template<typename _val>
	inline bool	is_low_freq(const _val *seq) const
	{
		double f = 0;
		for(unsigned i=0;i<weight_;++i) {
			_val l = seq[positions_[i]];
			if(l == Value_traits<_val>::MASK_CHAR || l == String_set<_val>::PADDING_CHAR)
				return false;
			l = mask_critical(l);
			unsigned r = Reduction<_val>::reduction(l);
			f += background_freq[r];
		}
		return !use_seed_freq<_val>() || f <= program_options::max_seed_freq;
	}

	template<typename _val>
	inline bool	is_low_freq_rev(const _val *seq) const
	{
		double f = 0;
		for(unsigned i=0;i<weight_;++i) {
			_val l = seq[(int)positions_[i]-(int)length_];
			if(l == Value_traits<_val>::MASK_CHAR || l == String_set<_val>::PADDING_CHAR)
				return false;
			l = mask_critical(l);
			unsigned r = Reduction<_val>::reduction(l);
			f += background_freq[r];
		}
		return !use_seed_freq<_val>() || f <= program_options::max_seed_freq;
	}

	uint32_t length_, weight_, positions_[Const::max_seed_weight], d_, mask_, rev_mask_, id_;

};

#endif /* SHAPE_H_ */

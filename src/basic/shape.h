/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen
						
Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#pragma once
#include <string.h>
#include "const.h"
#include "value.h"
#include "seed.h"
#include "../stats/score_matrix.h"
#include "reduction.h"
#include "config.h"

struct Shape
{

	Shape():
		length_ (0),
		weight_ (0),
		d_ (0),
		mask_ (0),
		rev_mask_ (0)
	{ memset(positions_, 0, sizeof(uint32_t)*Const::max_seed_weight); }

	Shape(const char* code, unsigned id) :
		weight_(0),
		mask_(0),
		rev_mask_(0),
		long_mask_(0)
		//long_mask_sse_(_mm_setzero_si128())
	{
		assert(strlen(code) <= 32);
		memset(positions_, 0, sizeof(uint32_t)*Const::max_seed_weight);
		const uint64_t b = Reduction::reduction.bit_size();
		unsigned i (0);
		for(;i<strlen(code);++i) {
			rev_mask_ <<= 1;
			long_mask_ <<= b;
			//long_mask_sse_ = _mm_slli_si128(long_mask_sse_, 1);
			if(code[i] == '1') {
				assert(weight_ < Const::max_seed_weight);
				positions_[weight_] = i;
				++weight_;
				mask_ |= 1 << i;
				rev_mask_ |= 1;
				long_mask_ |= (1 << b) - 1;
				//long_mask_sse_ = _mm_insert_epi8(long_mask_sse_, 0xff, 0);
			}
		}
		length_ = i;
		d_ = positions_[weight_/2-1];
	}

	inline bool set_seed(PackedSeed &s, const Letter *seq) const
	{
		s = 0;
#ifdef FREQUENCY_MASKING
		double f = 0;
#endif
		for (int i = 0; i < weight_; ++i) {
			Letter l = seq[positions_[i]];
#ifdef SEQ_MASK
			l &= LETTER_MASK;
#endif
			if (!is_amino_acid(l))
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

	inline bool set_seed_shifted(PackedSeed &s, const Letter *seq) const
	{
		s = 0;
		const uint64_t b = Reduction::reduction.bit_size();
		for (int i = 0; i < weight_; ++i) {
			Letter l = seq[positions_[i]];
#ifdef SEQ_MASK
			l &= LETTER_MASK;
#endif
			if (l == value_traits.mask_char || l == Sequence::DELIMITER || l == STOP_LETTER)
				return false;
			unsigned r = Reduction::reduction(l);
			s <<= b;
			s |= uint64_t(r);
		}
		return true;
	}

	inline bool set_seed_reduced(PackedSeed &s, const Letter *seq) const
	{
		s = 0;
		for (int i = 0; i < weight_; ++i) {
			Letter l = seq[positions_[i]];
#ifdef SEQ_MASK
			l &= LETTER_MASK;
#endif
			if (l == MASK_LETTER)
				return false;
			s *= Reduction::reduction.size();
			s += uint64_t(l);
		}
		return true;
	}

	inline bool set_seed(Seed &s, const Letter *seq) const
	{
		for (int i = 0; i < weight_; ++i) {
			Letter l = seq[positions_[i]];
#ifdef SEQ_MASK
			l &= LETTER_MASK;
#endif
			if (l >= 20)
				return false;
			s[i] = l;
		}
		return true;
	}
	
	friend std::ostream& operator<<(std::ostream&s, const Shape &sh)
	{
		for (int i = 0; i < sh.length_; ++i)
			s << ((sh.mask_ & (1 << i)) ? '1' : '0');
		return s;
	}

	bool contiguous() const
	{
		return length_ == weight_;
	}

	uint64_t long_mask() const
	{
		return long_mask_;
	}

	int32_t length_, weight_, positions_[Const::max_seed_weight];
	uint32_t d_, mask_, rev_mask_;
	uint64_t long_mask_;
	//__m128i long_mask_sse_;

};

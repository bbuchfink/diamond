/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once
#include <string.h>
#include "const.h"
#include "value.h"
#include "seed.h"
#include "reduction.h"
#include "util/math/integer.h"
#include "util/system.h"

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
		const uint64_t b = Reduction::get_reduction().bit_size();
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
			unsigned r = Reduction::get_reduction()(l);
#ifdef FREQUENCY_MASKING
			f += background_freq[r];
#endif
			s *= Reduction::get_reduction().size();
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
		const uint64_t b = Reduction::get_reduction().bit_size();
		for (int i = 0; i < weight_; ++i) {
			Letter l = seq[positions_[i]];
#ifdef SEQ_MASK
			l &= LETTER_MASK;
#endif
			if (l == value_traits.mask_char || l == Sequence::DELIMITER || l == STOP_LETTER)
				return false;
			unsigned r = Reduction::get_reduction()(l);
			s <<= b;
			s |= uint64_t(r);
		}
		return true;
	}
	
	template<int W, typename It>
	inline bool set_seed_reduced(PackedSeed& s, It seq) const
	{
		static_assert(W >= 2, "W must be >= 2");
		const PackedSeed size = Reduction::get_reduction().size();

		Letter letters[W];
		unsigned bad = 0;
#pragma GCC unroll 16
		for (int i = 0; i < W; ++i) {
			Letter l = seq[positions_[i]];
#ifdef SEQ_MASK
			l &= LETTER_MASK;
#endif
			bad |= unsigned(l == MASK_LETTER);
			letters[i] = l;
		}
		if (UNLIKELY(bad))
			return false;

		const PackedSeed s2 = size * size;
		PackedSeed E = letters[0];
		PackedSeed O = letters[1];

#pragma GCC unroll 16
		for (int i = 2; i + 1 < W; i += 2) {
			E = E * s2 + letters[i];
			O = O * s2 + letters[i + 1];
		}

		if (W % 2 == 0) {
			s = E * size + O;
		}
		else {
			E = E * s2 + letters[W - 1];
			s = E + O * size;
		}
		return true;
	}
	
	template<typename It>
	inline bool set_seed_reduced(PackedSeed& s, It seq) const
	{
		switch (weight_) {
		case 7:
			return set_seed_reduced<7>(s, seq);
		case 8:
			return set_seed_reduced<8>(s, seq);
		case 9:
			return set_seed_reduced<9>(s, seq);
		case 10:
			return set_seed_reduced<10>(s, seq);
		case 11:
			return set_seed_reduced<11>(s, seq);
		case 12:
			return set_seed_reduced<12>(s, seq);
		}
		UNREACHABLE;
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

	int bit_length() const {
		return ::bit_length(power((int64_t)Reduction::get_reduction().size(), (int64_t)weight_) - 1);
	}

	int32_t length_, weight_, positions_[Const::max_seed_weight];
	uint32_t d_, mask_, rev_mask_;
	uint64_t long_mask_;
	//__m128i long_mask_sse_;

};
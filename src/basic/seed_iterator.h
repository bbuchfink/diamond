/****
DIAMOND protein aligner
Copyright (C) 2016-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
						
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
#include "shape.h"
#include "sequence.h"
#include "../util/hash_function.h"
#include "../util/algo/MurmurHash3.h"

struct Seed_iterator
{
	Seed_iterator(vector<Letter> &seq, const Shape &sh):
		ptr_ (seq.data()),
		end_ (ptr_ + seq.size() - sh.length_ + 1)
	{}
	bool good() const
	{
		return ptr_ < end_;
	}
	bool get(uint64_t &seed, const Shape &sh)
	{
		return sh.set_seed_reduced(seed, ptr_++);
	}
private:
	const Letter *ptr_, *end_;
};

template<uint64_t _b>
struct Hashed_seed_iterator2
{
	Hashed_seed_iterator2(const Sequence &seq, const Shape &sh):
		ptr_(seq.data()),
		end_(ptr_ + seq.length()),
		last_(0)
	{
		for (uint64_t i = 0; (i < sh.length_ - 1) && ptr_ < end_; ++i) {
#ifdef SEQ_MASK
			last_ = (last_ << _b) | Reduction::reduction(letter_mask(*(ptr_++)));
#else
			last_ = (last_ << _b) | Reduction::reduction(*(ptr_++));
#endif
		}
	}
	bool good() const
	{
		return ptr_ < end_;
	}
	bool get(uint64_t &seed, uint64_t mask)
	{
		last_ <<= _b;
#ifdef SEQ_MASK
		const Letter l = *(ptr_++) & LETTER_MASK;
#else
		const Letter l = *(ptr_++);
#endif
		if (!is_amino_acid(l))
			return false;
		last_ |= Reduction::reduction(l);
		seed = murmur_hash()(last_ & mask);
		return true;
	}
private:
	const Letter *ptr_, *end_;
	uint64_t last_;
};

static inline uint64_t hash(__m128i a) {
	char in[16], seed[16], out[16];
	_mm_storeu_si128((__m128i*)in, a);
	std::fill(seed, seed + 16, 0);
	MurmurHash3_x64_128((void*)in, 16, seed, (void*)out);
	return *(uint64_t*)out;
	const __m128i prime32 = _mm_set1_epi32((int)0x9E3779B1U);
	__m128i const acc_vec = a;
	__m128i const shifted = _mm_srli_epi64(acc_vec, 47);
	__m128i const data_key = _mm_xor_si128(acc_vec, shifted);
	__m128i const data_key_hi = _mm_shuffle_epi32(data_key, _MM_SHUFFLE(0, 3, 0, 1));
	__m128i const prod_lo = _mm_mul_epu32(data_key, prime32);
	__m128i const prod_hi = _mm_mul_epu32(data_key_hi, prime32);
	return _mm_extract_epi64(_mm_add_epi64(prod_lo, _mm_slli_epi64(prod_hi, 32)), 0);
}

template<uint64_t _b>
struct Hashed_seed_iterator
{
	Hashed_seed_iterator(const Sequence& seq, const Shape& sh) :
		ptr_(seq.data()),
		end_(ptr_ + seq.length()),
		last_(_mm_setzero_si128())
	{
		for (uint64_t i = 0; (i < sh.length_ - 1) && ptr_ < end_; ++i) {
#ifdef SEQ_MASK
			last_ = _mm_slli_si128(last_, 1);
			const unsigned l = Reduction::reduction(letter_mask(*(ptr_++)));
			last_ = _mm_insert_epi8(last_, l, 0);
#else
			last_ = (last_ << _b) | Reduction::reduction(*(ptr_++));
#endif
		}
	}
	bool good() const
	{
		return ptr_ < end_;
	}
	bool get(uint64_t& seed, __m128i mask)
	{
		last_ = _mm_slli_si128(last_, 1);
#ifdef SEQ_MASK
		const Letter l = *(ptr_++) & LETTER_MASK;
#else
		const Letter l = *(ptr_++);
#endif
		if (!is_amino_acid(l))
			return false;
		last_ = _mm_insert_epi8(last_, Reduction::reduction(l), 0);
		seed = hash(_mm_and_si128(last_, mask));
		return true;
	}
private:
	const Letter* ptr_, * end_;
	__m128i last_;
};

template<uint64_t _l, uint64_t _b>
struct Contiguous_seed_iterator
{
	Contiguous_seed_iterator(const Sequence &seq) :
		ptr_(seq.data()),
		end_(ptr_ + seq.length()),
		last_(0)
	{
		for (uint64_t i = 0; i < _l - 1; ++i)
#ifdef SEQ_MASK
			last_ = (last_ << _b) | Reduction::reduction(letter_mask(*(ptr_++)));
#else
			last_ = (last_ << _b) | Reduction::reduction(*(ptr_++));
#endif
	}
	bool good() const
	{
		return ptr_ < end_;
	}
	bool get(uint64_t &seed)
	{
		last_ <<= _b;
		last_ &= (1 << (_b*_l)) - 1;
#ifdef SEQ_MASK
		const Letter l = *(ptr_++) & LETTER_MASK;
#else
		const Letter l = *(ptr_++);
#endif
		if (l == value_traits.mask_char)
			return false;
		last_ |= Reduction::reduction(l);
		seed = last_;
		return true;
	}
	static uint64_t length()
	{
		return _l;
	}
private:
	const Letter *ptr_, *end_;
	uint64_t last_;
};

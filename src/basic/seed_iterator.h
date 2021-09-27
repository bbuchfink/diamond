/****
DIAMOND protein aligner
Copyright (C) 2016-2021 Max Planck Society for the Advancement of Science e.V.
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

struct SeedIterator
{
	SeedIterator(vector<Letter> &seq, const Shape &sh):
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
	SeedIterator& operator++() {
		ptr_++;
		return *this;
	}
private:
	const Letter *ptr_, *end_;
};

template<uint64_t B>
struct HashedSeedIterator
{
	HashedSeedIterator(const Sequence &seq, const Shape &sh):
		ptr_(seq.data()),
		end_(ptr_ + seq.length()),
		last_(0)
	{
		for (int i = 0; (i < sh.length_ - 1) && ptr_ < end_; ++i) {
			last_ = (last_ << B) | Reduction::reduction(letter_mask(*(ptr_++)));
		}
	}
	bool good() const
	{
		return ptr_ < end_;
	}
	bool get(uint64_t &seed, uint64_t mask)
	{
		last_ <<= B;
		const Letter l = letter_mask(*(ptr_++));
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

struct FilterMaskedSeeds { };

template<int L, uint64_t B, typename Filter>
struct ContiguousSeedIterator
{
	ContiguousSeedIterator(const Sequence &seq) :
		ptr_(seq.data()),
		end_(ptr_ + seq.length()),
		last_(0)
	{
		for (int i = 0; i < L - 1; ++i)
			last_ = (last_ << B) | Reduction::reduction(letter_mask(*(ptr_++)));
	}
	bool good() const
	{
		return ptr_ < end_;
	}
	bool get(uint64_t &seed)
	{
		for (;;) {
			last_ <<= B;
			last_ &= (uint64_t(1) << (B*L)) - 1;
			const Letter l = letter_mask(*(ptr_++));
			last_ |= Reduction::reduction(l);
			seed = last_;
			return true;
		}
	}
	static int length()
	{
		return L;
	}
private:
	const Letter *ptr_, *end_;
	uint64_t last_;
	unsigned mask_;
};

template<int L, uint64_t B>
struct ContiguousSeedIterator<L, B, FilterMaskedSeeds>
{
	ContiguousSeedIterator(const Sequence &seq) :
		ptr_(seq.data()),
		end_(ptr_ + seq.length()),
		last_(0),
		mask_(0)
	{
		for (int i = 0; i < L - 1; ++i) {
			const Letter l = letter_mask(*(ptr_++));
			last_ = (last_ << B) | Reduction::reduction(l);
			if (!is_amino_acid(l))
				mask_ |= 1;
			mask_ <<= 1;
		}
	}
	bool good() const
	{
		return ptr_ < end_;
	}
	bool get(uint64_t &seed)
	{
		for (;;) {
			last_ <<= B;
			last_ &= (uint64_t(1) << (B*L)) - 1;
			mask_ <<= 1;
			mask_ &= (1 << L) - 1;
			const Letter l = letter_mask(*(ptr_++));
			const unsigned r = Reduction::reduction(l);
			last_ |= r;
			seed = last_;
			if (!is_amino_acid(l))
				mask_ |= 1;
			return mask_ == 0;
		}
	}
	static int length()
	{
		return L;
	}
private:
	const Letter *ptr_, *end_;
	uint64_t last_;
	unsigned mask_;
};
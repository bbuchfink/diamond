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

#ifndef SEED_ITERATOR_H_
#define SEED_ITERATOR_H_

#include "shape.h"
#include "sequence.h"
#include "../util/hash_function.h"

struct Seed_iterator
{
	Seed_iterator(vector<char> &seq, const shape &sh):
		ptr_ (seq.data()),
		end_ (ptr_ + seq.size() - sh.length_ + 1)
	{}
	bool good() const
	{
		return ptr_ < end_;
	}
	bool get(uint64_t &seed, const shape &sh)
	{
		return sh.set_seed_reduced(seed, ptr_++);
	}
private:
	const char *ptr_, *end_;
};

template<uint64_t _b>
struct Hashed_seed_iterator
{
	Hashed_seed_iterator(const sequence &seq, const shape &sh):
		ptr_(seq.data()),
		end_(ptr_ + seq.length()),
		last_(0)
	{
		for (uint64_t i = 0; i < sh.length_ - 1; ++i)
			last_ = (last_ << _b) | Reduction::reduction(*(ptr_++));
	}
	bool good() const
	{
		return ptr_ < end_;
	}
	bool get(uint64_t &seed, uint64_t mask)
	{
		last_ <<= _b;
		const char l = *(ptr_++);
		if (l == value_traits.mask_char)
			return false;
		last_ |= Reduction::reduction(l);
		seed = murmur_hash()(last_ & mask);
		return true;
	}
private:
	const char *ptr_, *end_;
	uint64_t last_;
};

template<uint64_t _l, uint64_t _b>
struct Contiguous_seed_iterator
{
	Contiguous_seed_iterator(const sequence &seq) :
		ptr_(seq.data()),
		end_(ptr_ + seq.length()),
		last_(0)
	{
		for (uint64_t i = 0; i < _l - 1; ++i)
			last_ = (last_ << _b) | Reduction::reduction(*(ptr_++));
	}
	bool good() const
	{
		return ptr_ < end_;
	}
	bool get(uint64_t &seed)
	{
		last_ <<= _b;
		last_ &= (1 << (_b*_l)) - 1;
		const char l = *(ptr_++);
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
	const char *ptr_, *end_;
	uint64_t last_;
};

#endif
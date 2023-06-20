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
#include <deque>
#include "shape.h"
#include "sequence.h"
#include "../util/hash_function.h"
#include "../util/algo/MurmurHash3.h"

struct SeedIterator
{
	SeedIterator(std::vector<Letter> &seq, const Shape &sh):
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

struct MinimizerIterator
{
	MinimizerIterator(std::vector<Letter>& seq, const Shape& sh, Loc window) :
		ptr_(seq.data()),
		begin_(seq.data()),
		end_(ptr_ + seq.size() - sh.length_ + 1),
		window_(window),
		sh_(sh)
	{
		next();
		if (good())
			min_idx_ = get();
	}
	bool good() const
	{
		return (Loc)seeds_.size() == window_;
	}
	uint64_t operator*() const {
		return seeds_[min_idx_];
	}
	MinimizerIterator& operator++() {
		int m = 0;
		const uint64_t current = **this;
		do {
			seeds_.pop_front();
			hashes_.pop_front();
			pos_.pop_front();
			next();
		} while (good() && seeds_[m = get()] == current);
		min_idx_ = m;
		return *this;
	}
	Loc pos() const {
		return pos_[min_idx_];
	}
private:
	void next() {
		while ((Loc)seeds_.size() < window_ && ptr_ < end_) {
			uint64_t s;
			if (sh_.set_seed_reduced(s, ptr_)) {
				seeds_.push_back(s);
				hashes_.push_back(MurmurHash()(s));
				pos_.push_back(Loc(ptr_ - begin_));
			}
			++ptr_;
		}
	}
	int get() const {
		std::deque<uint64_t>::const_iterator i = hashes_.begin(), j = i;
		uint64_t s = *i++;
		while (i < hashes_.end()) {
			if (*i < s) {
				s = *i;
				j = i;
			}
			++i;
		}
		return int(j - hashes_.begin());
	}
	const Letter* ptr_, *begin_, *end_;
	std::deque<uint64_t> seeds_, hashes_;
	std::deque<Loc> pos_;
	const Loc window_;
	const Shape& sh_;
	int min_idx_;
};

struct SketchIterator
{
	SketchIterator(std::vector<Letter>& seq, const Shape& sh, Loc n)
	{
		std::vector<Kmer> v;
		v.reserve(seq.size() - sh.length_ + 1);
		const Letter* end = seq.data() + seq.size() - sh.length_ + 1;
		uint64_t s;
		for (const Letter* p = seq.data(); p < end; ++p)
			if (sh.set_seed_reduced(s, p))
				v.emplace_back(s, MurmurHash()(s), Loc(p - seq.data()));
		std::sort(v.begin(), v.end());
		data_.insert(data_.end(), v.begin(), v.begin() + std::min(n, (Loc)v.size()));
		it_ = data_.begin();
	}
	bool good() const {
		return it_ < data_.end();
	}
	uint64_t operator*() const {
		return it_->seed;
	}
	Loc pos() const {
		return it_->pos;
	}
	SketchIterator& operator++() {
		++it_;
		return *this;
	}
private:
	struct Kmer {
		Kmer(uint64_t seed, uint64_t hash, Loc pos) : seed(seed), hash(hash), pos(pos) {}
		uint64_t seed, hash;
		Loc pos;
		bool operator<(const Kmer& k) const {
			return hash < k.hash;
		}
	};
	std::vector<Kmer> data_;
	std::vector<Kmer>::const_iterator it_;
};

template<uint64_t B>
struct HashedSeedIterator
{
	HashedSeedIterator(Letter* seq, Loc len, const Shape &sh):
		long_mask(sh.long_mask()),
		ptr_(seq),
		end_(ptr_ + len),
		last_(0)
	{
		for (int i = 0; i < sh.length_ && ptr_ < end_; ++i)
			last_ = (last_ << B) | Reduction::reduction(letter_mask(*(ptr_++)));
	}
	bool good() const
	{
		return ptr_ <= end_;
	}
	uint64_t operator*() const {
		return MurmurHash()(last_ & long_mask);
	}
	HashedSeedIterator& operator++() {
		while (ptr_ < end_) {
			last_ <<= B;
			const Letter l = letter_mask(*(ptr_++));
			if (!is_amino_acid(l))
				continue;
			last_ |= Reduction::reduction(l);
			return *this;
		}
		++ptr_;
	}
	Letter* seq_ptr(const Shape& sh) const {
		return ptr_ - sh.length_;
	}
private:
	const uint64_t long_mask;
	Letter *ptr_, *end_;
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
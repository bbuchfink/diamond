/****
DIAMOND protein aligner
Copyright (C) 2022 Max Planck Society for the Advancement of Science e.V.

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
#include <stdint.h>
#include "../../basic/value.h"
#include "../../basic/sequence.h"
#include "../math/integer.h"

struct IdentityReduction {
	int bit_size() const {
		return 5;
	}
	uint64_t size() const {
		return 20;
	}
	uint64_t operator()(uint64_t x) const {
		return x;
	}
};

template<size_t K>
struct Kmer {
	Kmer() :
		code(0)
	{ }
	Kmer(const char* s) :
		code(0)
	{
		assert(strlen(s) == K);
		for (size_t i = 0; i < K; ++i)
			code = (code * TRUE_AA) + (uint64_t)amino_acid_traits.from_char(*s++);
	}
	operator uint64_t() const {
		return code;
	}
	struct Hash {
		size_t operator()(const Kmer<K> k) const {
			return std::hash<uint64_t>()(k.code);
		}
	};
	uint64_t code;
};

template<size_t K, typename R = IdentityReduction>
struct KmerIterator {
	KmerIterator(const Sequence& seq, const R& reduction = R()) :
		reduction_(reduction),
		ptr_(seq.data() - 1),
		end_(seq.end()),
		mod_(power((uint64_t)reduction_.size(), K - 1))
	{
		inc(0, 1);
	}
	Kmer<K> operator*() const {
		return kmer_;
	}
	bool good() const {
		return ptr_ < end_;
	}
	KmerIterator& operator++() {
		inc(K - 1, mod_);
		return *this;
	}
	Loc operator-(const Letter* ptr) {
		return Loc(ptr_ + 1 - K - ptr);
	}
private:
	void inc(uint64_t n, uint64_t mod) {
		kmer_.code %= mod;
		do {
			++ptr_;
			if (ptr_ == end_)
				return;
			const uint64_t l = letter_mask(*ptr_);
			if (l < TRUE_AA) {
				kmer_.code = (kmer_.code * reduction_.size()) + reduction_(l);
				++n;
			}
			else {
				kmer_.code = 0;
				n = 0;
			}
		} while (n < K);
		//kmer_.code &= (uint64_t(1) << (reduction_.size() * K)) - 1;
	}
	const R reduction_;
	const Letter* ptr_, *const end_;
	const uint64_t mod_;
	Kmer<K> kmer_;
};

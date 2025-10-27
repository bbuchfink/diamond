/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <stdint.h>
#include "basic/value.h"
#include "basic/sequence.h"
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
		mod_(power((size_t)reduction_.size(), K - 1))
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

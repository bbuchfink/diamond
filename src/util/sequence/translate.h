/****
Copyright © 2013-2026 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-2-Clause

#pragma once
#include <vector>
#include "basic/sequence.h"
#include "basic/value.h"

struct Translator
{

public:

	static const Letter reverseLetter[5];
	static Letter lookup[5][5][5];
	static Letter lookupReverse[5][5][5];
	static const char* codes[27];

	static void init(unsigned id);

	static Letter getReverseComplement(Letter letter)
	{ return reverseLetter[(int) letter]; }

	static Letter getAminoAcid(const Sequence& dnaSequence, size_t pos)
	{ return lookup[(int) dnaSequence[pos]][(int)dnaSequence[pos+1]][(int)dnaSequence[pos+2]]; }

	static Letter getAminoAcidReverse(const Sequence&  dnaSequence, size_t pos)
	{ return lookupReverse[(int) dnaSequence[pos + 2]][(int)dnaSequence[pos + 1]][(int)dnaSequence[pos]]; }

	static std::vector<Letter> reverse(const Sequence &seq)
	{
		std::vector<Letter> r;
        r.reserve(seq.len_+1);
		for(const Letter* end = seq.data_ + seq.len_ -1 ; end >= seq.data_; --end)
			r.push_back(getReverseComplement(*end));
		return r;
	}

	static size_t translate(const Sequence& dnaSequence, std::array<std::vector<Letter>, 6>& proteins)
	{
		const size_t length_ = dnaSequence.length();
		if (length_ < 3) {
			for (int i = 0; i < 6; ++i)
				proteins[i].clear();
			return 0;
		}
		const size_t l1 = length_ / 3;
		proteins[0].resize(l1);
		proteins[3].resize(l1);
		const size_t l2 = (length_ - 1) / 3;
		proteins[1].resize(l2);
		proteins[4].resize(l2);
		const size_t l3 = (length_ - 2) / 3;
		proteins[2].resize(l3);
		proteins[5].resize(l3);

		size_t r = length_ - 2;
		unsigned pos = 0;
		unsigned i = 0;
		while(r > 2) {
			proteins[0][i] = getAminoAcid(dnaSequence, pos++);
			proteins[3][i] = getAminoAcidReverse(dnaSequence, --r);
			proteins[1][i] = getAminoAcid(dnaSequence, pos++);
			proteins[4][i] = getAminoAcidReverse(dnaSequence, --r);
			proteins[2][i] = getAminoAcid(dnaSequence, pos++);
			proteins[5][i] = getAminoAcidReverse(dnaSequence, --r);
			++i;
		}
		if(r) {
			proteins[0][i] = getAminoAcid(dnaSequence, pos++);
			proteins[3][i] = getAminoAcidReverse(dnaSequence, --r);
		}
		if(r) {
			proteins[1][i] = getAminoAcid(dnaSequence, pos);
			proteins[4][i] = getAminoAcidReverse(dnaSequence, r - 1);
		}
		return 2 * (l1 + l2 + l3);
	}

	static Letter const* nextChar(Letter const*p, Letter const*end)
	{
		while(*(p) != STOP_LETTER && p < end)
			++p;
		return p;
	}

	static void mask_runs(std::vector<Letter> &query, unsigned run_len)
	{
		Letter *last = &query[0]-1, *i = &query[0], *end = &query.back();
		while (i <= end) {
			if(*i == STOP_LETTER) {
				if(last != 0 && i - last - 1 < (int)run_len) {
					for(Letter *j = last+1; j < i; ++j)
						*j = 23;
				}
				last = i;
			}
			++i;
		}
		if(last != 0 && i - last - 1 < (int)run_len) {
			for(Letter *j = last+1; j < i; ++j)
				*j = 23;
		}
	}

	static unsigned computeGoodFrames(std::vector<Letter>  const *queries, unsigned runLen)
	{
		unsigned set = 0;

		for (unsigned i = 0; i < 6; ++i) {
			if (queries[i].size() > 0) {
				unsigned run = 0;
				Letter const*p =  &(queries[i][0]);
				Letter const*q;
				Letter const*end = p + queries[i].size();
				while((q = nextChar(p, end)) < end) {
					run = (unsigned)(q-p);
					if (run >= runLen)
						set |= 1 << i;
					p=q+1;
				}
				run = (unsigned)(q-p);
				if (run >= runLen)
					set |= 1 << i;
			}
		}
		return set;
	}

	static void mask_runs(std::vector<Letter> *queries, unsigned run_len)
	{
		for (unsigned i = 0; i < 6; ++i)
			mask_runs(queries[i], run_len);
	}

};


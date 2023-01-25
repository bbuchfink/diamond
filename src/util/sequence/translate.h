/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

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
#include <vector>
#include "../../basic/value.h"

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

	static std::vector<Letter> reverse(const std::vector<Letter> &seq)
	{
		std::vector<Letter> r;
		for(std::vector<Letter>::const_reverse_iterator i=seq.rbegin(); i!=seq.rend(); ++i)
			r.push_back(getReverseComplement(*i));
		return r;
	}

	static size_t translate(const Sequence& dnaSequence, std::vector<Letter> *proteins)
	{
		size_t length_ = dnaSequence.length(), d, n;
		proteins[0].resize(d = length_ / 3);
		proteins[3].resize(d);
		n = 2*d;
		proteins[1].resize(d = (length_-1) / 3);
		proteins[4].resize(d);
		n += 2*d;
		proteins[2].resize(d = (length_-2) / 3);
		proteins[5].resize(d);
		n += 2*d;

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
			proteins[4][i] = getAminoAcidReverse(dnaSequence, r);
		}
		return n;
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


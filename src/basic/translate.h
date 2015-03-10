/****
Copyright (c) 2014, University of Tuebingen
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****
Author: Benjamin Buchfink
****/

#ifndef TRANSLATE_H_
#define TRANSLATE_H_

struct Translator
{

public:

	static const Nucleotide reverseNucleotide[5];
	static const Amino_acid lookup[5][5][5];
	static const Amino_acid lookupReverse[5][5][5];
	static const Amino_acid STOP;

	static Nucleotide getReverseComplement(Nucleotide nucleotide)
	{ return reverseNucleotide[(int) nucleotide]; }

	static Amino_acid getAminoAcid(vector<Nucleotide> const &dnaSequence, size_t pos)
	{ return lookup[(int) dnaSequence[pos]][(int)dnaSequence[pos+1]][(int)dnaSequence[pos+2]]; }

	static Amino_acid getAminoAcidReverse(vector<Nucleotide> const &dnaSequence, size_t pos)
	{ return lookupReverse[(int) dnaSequence[pos + 2]][(int)dnaSequence[pos + 1]][(int)dnaSequence[pos]]; }

	static vector<Nucleotide> reverse(const vector<Nucleotide> &seq)
	{
		vector<Nucleotide> r;
		for(vector<Nucleotide>::const_reverse_iterator i=seq.rbegin(); i!=seq.rend(); ++i)
			r.push_back(getReverseComplement(*i));
		return r;
	}

	static size_t translate(vector<Nucleotide> const &dnaSequence, vector<Amino_acid> *proteins)
	{
		size_t length_ = dnaSequence.size(), d, n;
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

	static Amino_acid const* nextChar(Amino_acid const*p, Amino_acid const*end)
	{
		while(*(p) != STOP && p < end)
			++p;
		return p;
	}

	static void mask_runs(vector<Amino_acid> &query, unsigned run_len)
	{
		Amino_acid *last = &query[0]-1, *i = &query[0], *end = &query.back();
		while (i <= end) {
			if(*i == STOP) {
				if(last != 0 && i - last - 1 < run_len) {
					for(Amino_acid *j = last+1; j < i; ++j)
						*j = 23;
				}
				last = i;
			}
			++i;
		}
		if(last != 0 && i - last - 1 < run_len) {
			for(Amino_acid *j = last+1; j < i; ++j)
				*j = 23;
		}
	}

	static unsigned computeGoodFrames(vector<Amino_acid>  const *queries, unsigned runLen)
	{
		unsigned set = 0;

		for (unsigned i = 0; i < 6; ++i) {
			if (queries[i].size() > 0) {
				unsigned run = 0;
				Amino_acid const*p =  &(queries[i][0]);
				Amino_acid const*q;
				Amino_acid const*end = p + queries[i].size();
				while((q = nextChar(p, end)) < end) {
					run = q-p;
					if (run >= runLen)
						set |= 1 << i;
					p=q+1;
				}
				run = q-p;
				if (run >= runLen)
					set |= 1 << i;
			}
		}
		return set;
	}

	static void mask_runs(vector<Amino_acid> *queries, unsigned run_len)
	{
		for (unsigned i = 0; i < 6; ++i)
			mask_runs(queries[i], run_len);
	}

};

const Nucleotide Translator::reverseNucleotide[5] = { 3, 2, 1, 0, 4 };

const Amino_acid Translator::lookup[5][5][5] = {
{ { 11,2,11,2,23 },
{ 16,16,16,16,16 },
{ 1,15,1,15,23 },
{ 9,9,12,9,23 },
{ 23,23,23,23,23 },
 },
{ { 5,8,5,8,23 },
{ 14,14,14,14,14 },
{ 1,1,1,1,1 },
{ 10,10,10,10,10 },
{ 23,23,23,23,23 },
 },
{ { 6,3,6,3,23 },
{ 0,0,0,0,0 },
{ 7,7,7,7,7 },
{ 19,19,19,19,19 },
{ 23,23,23,23,23 },
 },
{ { 23,18,23,18,23 },
{ 15,15,15,15,15 },
{ 23,4,17,4,23 },
{ 10,13,10,13,23 },
{ 23,23,23,23,23 },
 },
{ { 23,23,23,23,23 },
{ 23,23,23,23,23 },
{ 23,23,23,23,23 },
{ 23,23,23,23,23 },
{ 23,23,23,23,23 },
} };

const Amino_acid Translator::lookupReverse[5][5][5] = {
{ { 13,10,13,10,23 },
{ 4,17,4,23,23 },
{ 15,15,15,15,15 },
{ 18,23,18,23,23 },
{ 23,23,23,23,23 },
 },
{ { 19,19,19,19,19 },
{ 7,7,7,7,7 },
{ 0,0,0,0,0 },
{ 3,6,3,6,23 },
{ 23,23,23,23,23 },
 },
{ { 10,10,10,10,10 },
{ 1,1,1,1,1 },
{ 14,14,14,14,14 },
{ 8,5,8,5,23 },
{ 23,23,23,23,23 },
 },
{ { 9,12,9,9,23 },
{ 15,1,15,1,23 },
{ 16,16,16,16,16 },
{ 2,11,2,11,23 },
{ 23,23,23,23,23 },
 },
{ { 23,23,23,23,23 },
{ 23,23,23,23,23 },
{ 23,23,23,23,23 },
{ 23,23,23,23,23 },
{ 23,23,23,23,23 },
}};

const Amino_acid Translator::STOP (Value_traits<Amino_acid>::from_char('*'));

#endif /* TRANSLATE_H_ */

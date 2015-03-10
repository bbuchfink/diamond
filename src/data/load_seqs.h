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

#ifndef LOAD_SEQS_H_
#define LOAD_SEQS_H_

#include "sequence_set.h"
#include "../basic/translate.h"
#include "../util/seq_file_format.h"

struct Single_strand { };
struct Double_strand { };

template<typename _ival, typename _val, typename _strand>
size_t push_seq(String_set<_val> &ss, String_set<Nucleotide>& source_seqs, const vector<_ival> &seq)
{ ss.push_back(seq); return seq.size(); }

template<>
size_t push_seq<Amino_acid,Nucleotide,Single_strand>(String_set<Nucleotide> &ss, String_set<Nucleotide>& source_seqs, const vector<Amino_acid> &seq)
{ return 0; }

template<>
size_t push_seq<Nucleotide,Amino_acid,Double_strand>(String_set<Amino_acid> &ss, String_set<Nucleotide>& source_seqs, const vector<Nucleotide> &seq)
{
	source_seqs.push_back(seq);
	if(seq.size() < 2) {
		for(unsigned j=0;j<6;++j)
			ss.fill(0, Value_traits<Amino_acid>::MASK_CHAR);
		return 0;
	}
	vector<Amino_acid> proteins[6];
	size_t n = Translator::translate(seq, proteins);

	unsigned bestFrames (Translator::computeGoodFrames(proteins, program_options::get_run_len(seq.size()/3)));
	for(unsigned j = 0; j < 6; ++j) {
		if(bestFrames & (1 << j))
			ss.push_back(proteins[j]);
		else
			ss.fill(proteins[j].size(), Value_traits<Amino_acid>::MASK_CHAR);
	}
	return n;
}

template<>
size_t push_seq<Nucleotide,Nucleotide,Double_strand>(String_set<Nucleotide> &ss, String_set<Nucleotide>& source_seqs, const vector<Nucleotide> &seq)
{
	ss.push_back(seq);
	ss.push_back(Translator::reverse(seq));
	return seq.size()*2;
}

template<typename _ival, typename _val, typename _strand>
size_t load_seqs(Input_stream &file,
		const Sequence_file_format<_ival> &format,
		Sequence_set<_val>*& seqs,
		String_set<char,0>*& ids,
		Sequence_set<Nucleotide>*& source_seqs,
		size_t max_letters)
{
	seqs = new Sequence_set<_val> ();
	ids = new String_set<char,0> ();
	source_seqs = new Sequence_set<Nucleotide> ();
	size_t letters = 0, n = 0;
	vector<_ival> seq;
	vector<char> id;
	while(letters < max_letters && format.get_seq(id, seq, file)) {
		ids->push_back(id);
		letters += push_seq<_ival,_val,_strand>(*seqs, *source_seqs, seq);
		++n;
	}
	ids->finish_reserve();
	seqs->finish_reserve();
	source_seqs->finish_reserve();
	if(n == 0) {
		delete seqs;
		delete ids;
		delete source_seqs;
	}
	return n;
}

#endif /* LOAD_SEQS_H_ */

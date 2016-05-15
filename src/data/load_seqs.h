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

#include <iostream>
#include "sequence_set.h"
#include "../basic/translate.h"
#include "../util/seq_file_format.h"

size_t push_seq(Sequence_set &ss, Sequence_set& source_seqs, const vector<Letter> &seq)
{
	if (config.command == Config::blastp || config.command == Config::makedb) {
		ss.push_back(seq);
		return seq.size();
	}
	else {
		source_seqs.push_back(seq);
		if (seq.size() < 2) {
			for (unsigned j = 0; j<6; ++j)
				ss.fill(0, value_traits.mask_char);
			return 0;
		}
		vector<Letter> proteins[6];
		size_t n = Translator::translate(seq, proteins);

		unsigned bestFrames(Translator::computeGoodFrames(proteins, config.get_run_len((unsigned)seq.size() / 3)));
		for (unsigned j = 0; j < 6; ++j) {
			if (bestFrames & (1 << j))
				ss.push_back(proteins[j]);
			else
				ss.fill(proteins[j].size(), value_traits.mask_char);
		}
		return n;
	}
}

size_t load_seqs(Compressed_istream &file,
		const Sequence_file_format &format,
		Sequence_set** seqs,
		String_set<0>*& ids,
		Sequence_set*& source_seqs,
		size_t max_letters)
{
	*seqs = new Sequence_set ();
	ids = new String_set<0> ();
	source_seqs = new Sequence_set ();
	size_t letters = 0, n = 0;
	vector<Letter> seq;
	vector<char> id;
	try {
		while(letters < max_letters && format.get_seq(id, seq, file)) {
			ids->push_back(id);
			letters += push_seq(**seqs, *source_seqs, seq);
			++n;
		}
	} catch(invalid_sequence_char_exception &e) {
		std::cerr << n << endl;
		throw e;
	}
	ids->finish_reserve();
	(*seqs)->finish_reserve();
	source_seqs->finish_reserve();
	if(n == 0) {
		delete *seqs;
		delete ids;
		delete source_seqs;
	}
	return n;
}

#endif /* LOAD_SEQS_H_ */

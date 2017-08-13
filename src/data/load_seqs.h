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

#ifndef LOAD_SEQS_H_
#define LOAD_SEQS_H_

#include <iostream>
#include "sequence_set.h"
#include "../basic/translate.h"
#include "../util/seq_file_format.h"

inline size_t push_seq(Sequence_set &ss, Sequence_set** source_seqs, const vector<Letter> &seq, unsigned frame_mask)
{
	if (config.command == Config::blastp || config.command == Config::makedb || config.command == Config::random_seqs) {
		ss.push_back(seq);
		return seq.size();
	}
	else {
		(*source_seqs)->push_back(seq);
		if (seq.size() < 2) {
			for (unsigned j = 0; j<6; ++j)
				ss.fill(0, value_traits.mask_char);
			return 0;
		}
		vector<Letter> proteins[6];
		size_t n = Translator::translate(seq, proteins);

		unsigned bestFrames(Translator::computeGoodFrames(proteins, config.get_run_len((unsigned)seq.size() / 3)));
		for (unsigned j = 0; j < 6; ++j) {
			if ((bestFrames & (1 << j)) && (frame_mask & (1 << j)))
				ss.push_back(proteins[j]);
			else
				ss.fill(proteins[j].size(), value_traits.mask_char);
		}
		return n;
	}
}

inline size_t load_seqs(Input_stream &file,
	const Sequence_file_format &format,
	Sequence_set** seqs,
	String_set<0>*& ids,
	Sequence_set** source_seqs,
	size_t max_letters,
	const string &filter)
{
	*seqs = new Sequence_set();
	ids = new String_set<0>();
	if(source_seqs)
		*source_seqs = new Sequence_set();
	size_t letters = 0, n = 0;
	vector<Letter> seq;
	vector<char> id;
	string id2;

	unsigned frame_mask = (1 << 6) - 1;
	if (config.query_strands == "plus")
		frame_mask = (1 << 3) - 1;
	else if (config.query_strands == "minus")
		frame_mask = ((1 << 3) - 1) << 3;

	while (letters < max_letters && format.get_seq(id, seq, file)) {
		if (seq.size() > 0 && (filter.empty() || id2.assign(id.data(), id.data() + id.size()).find(filter, 0) != string::npos)) {
			ids->push_back(id);
			letters += push_seq(**seqs, source_seqs, seq, frame_mask);
			++n;
			if ((*seqs)->get_length() >(size_t)std::numeric_limits<int>::max())
				throw std::runtime_error("Number of sequences in file exceeds supported maximum.");
		}
	}
	ids->finish_reserve();
	(*seqs)->finish_reserve();
	if(source_seqs)
		(*source_seqs)->finish_reserve();
	if (n == 0) {
		delete *seqs;
		delete ids;
		if(source_seqs)
			delete *source_seqs;
	}
	return n;
}

#endif /* LOAD_SEQS_H_ */

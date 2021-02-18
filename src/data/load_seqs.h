/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen
						
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
#include <list>
#include "sequence_set.h"
#include "../basic/translate.h"
#include "../util/seq_file_format.h"

inline size_t push_seq(SequenceSet &ss, SequenceSet** source_seqs, const vector<Letter> &seq, unsigned frame_mask, Sequence_type seq_type)
{
	if (seq_type == Sequence_type::amino_acid) {
		ss.push_back(seq.cbegin(), seq.cend());
		return seq.size();
	}
	else {
		(*source_seqs)->push_back(seq.cbegin(), seq.cend());
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
				ss.push_back(proteins[j].cbegin(), proteins[j].cend());
			else
				ss.fill(proteins[j].size(), value_traits.mask_char);
		}
		return n;
	}
}

inline size_t load_seqs(std::list<TextInputFile>::iterator file_begin,
	std::list<TextInputFile>::iterator file_end,
	const Sequence_file_format &format,
	SequenceSet** seqs,
	String_set<char, '\0'>*& ids,
	SequenceSet** source_seqs,
	String_set<char, '\0'>** quals,
	size_t max_letters,
	const string &filter,
	const Value_traits &value_traits,
	size_t modulo = 1)
{
	*seqs = new SequenceSet();
	ids = new String_set<char, '\0'>();
	if(source_seqs)
		*source_seqs = new SequenceSet();
	if (quals)
		*quals = new String_set<char, '\0'>();
	size_t letters = 0, n = 0;
	vector<Letter> seq;
	string id;
	vector<char> qual;
	string id2;

	unsigned frame_mask = (1 << 6) - 1;
	if (config.query_strands == "plus")
		frame_mask = (1 << 3) - 1;
	else if (config.query_strands == "minus")
		frame_mask = ((1 << 3) - 1) << 3;

	std::list<TextInputFile>::iterator file_it = file_begin;
	bool read_success = true;

	while ((letters < max_letters || (n % modulo != 0)) && (read_success = format.get_seq(id, seq, *file_it, value_traits, quals ? &qual : nullptr))) {
		if (seq.size() > 0 && (filter.empty() || id2.assign(id.data(), id.data() + id.size()).find(filter, 0) != string::npos)) {
			ids->push_back(id.begin(), id.end());
			letters += push_seq(**seqs, source_seqs, seq, frame_mask, value_traits.seq_type);
			if (quals)
				(*quals)->push_back(qual.begin(), qual.end());
			++n;
			if ((*seqs)->get_length() > (size_t)std::numeric_limits<int>::max())
				throw std::runtime_error("Number of sequences in file exceeds supported maximum.");
		}
		++file_it;
		if (file_it == file_end)
			file_it = file_begin;
	}
	ids->finish_reserve();
	if (quals)
		(*quals)->finish_reserve();
	(*seqs)->finish_reserve();
	if(source_seqs)
		(*source_seqs)->finish_reserve();
	if (n == 0) {
		delete *seqs;
		delete ids;
		if(source_seqs)
			delete *source_seqs;
		if (quals)
			delete *quals;
	}
	if (file_it != file_begin || (!read_success && ++file_it != file_end && format.get_seq(id, seq, *file_it, value_traits, nullptr)))
		throw std::runtime_error("Unequal number of sequences in paired read files.");
	return n;
}

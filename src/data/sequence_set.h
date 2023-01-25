/****
DIAMOND protein aligner
Copyright (C) 2013-2021 Max Planck Society for the Advancement of Science e.V.
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
#include <string>
#include "../basic/sequence.h"
#include "string_set.h"
#include "../basic/value.h"

struct SequenceSet : public StringSetBase<Letter, Sequence::DELIMITER, 1>
{

	SequenceSet(Alphabet alphabet = Alphabet::STD);
	SequenceSet(StringSetBase&& string_set):
		StringSetBase(string_set),
		alphabet_(Alphabet::STD)
	{}
	
	void print_stats() const;

	Sequence operator[](size_t i) const
	{
		return Sequence(ptr(i), (Loc)length(i));
	}

	std::pair<Length, Length> len_bounds(Length min_len) const;

	Length max_len(size_t begin, size_t end) const;

	std::vector<size_t> partition(unsigned n_part) const;

	size_t reverse_translated_len(size_t i) const;

	TranslatedSequence translated_seq(const Sequence& source, size_t i) const;

	size_t avg_len() const;

	virtual ~SequenceSet();

	Alphabet alphabet() const {
		return alphabet_;
	}

	Alphabet& alphabet() {
		return alphabet_;
	}

	void convert_to_std_alph(size_t id);
	void convert_all_to_std_alph(size_t threads);
	std::vector<std::pair<Loc, BlockId>> lengths() const;
	
private:

	Alphabet alphabet_;

};

size_t max_id_len(const StringSet& ids);
std::vector<std::string> seq_titles(const char* title);
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
#include <string>
#include <algorithm>
#include <queue>
#include <thread>
#include "../basic/sequence.h"
#include "string_set.h"
#include "../basic/shape_config.h"
#include "../basic/seed_iterator.h"
#include "../util/ptr_vector.h"
#include "../basic/value.h"

struct Sequence_set : public String_set<Letter, Sequence::DELIMITER, 1>
{

	Sequence_set()
	{ }
	
	void print_stats() const
	{
		verbose_stream << "Sequences = " << this->get_length() << ", letters = " << this->letters() << ", average length = " << this->avg_len() << std::endl;
	}

	Sequence operator[](size_t i) const
	{
		return Sequence(ptr(i), length(i));
	}

	std::pair<size_t, size_t> len_bounds(size_t min_len) const
	{
		const size_t l(this->get_length());
		size_t max = 0, min = std::numeric_limits<size_t>::max();
		for (size_t i = 0; i < l; ++i) {
			max = std::max(this->length(i), max);
			min = this->length(i) >= min_len ? std::min(this->length(i), min) : min;
		}
		return std::pair<size_t, size_t>(min, max);
	}

	size_t max_len(size_t begin, size_t end) const
	{
		size_t max = 0;
		for (size_t i = begin; i < end; ++i)
			max = std::max(this->length(i), max);
		return max;
	}

	std::vector<size_t> partition(unsigned n_part) const
	{
		std::vector<size_t> v;
		const size_t l = (this->letters() + n_part - 1) / n_part;
		v.push_back(0);
		for (unsigned i = 0; i < this->get_length();) {
			size_t n = 0;
			while (i < this->get_length() && n < l)
				n += this->length(i++);
			v.push_back(i);
		}
		for (size_t i = v.size(); i < n_part + 1; ++i)
			v.push_back(this->get_length());
		return v;
	}

	size_t reverse_translated_len(size_t i) const
	{
		const size_t j(i - i % 6);
		const size_t l(this->length(j));
		if (this->length(j + 2) == l)
			return l * 3 + 2;
		else if (this->length(j + 1) == l)
			return l * 3 + 1;
		else
			return l * 3;
	}

	TranslatedSequence translated_seq(const Sequence &source, size_t i) const
	{
		if (!align_mode.query_translated)
			return TranslatedSequence((*this)[i]);
		return TranslatedSequence(source, (*this)[i], (*this)[i + 1], (*this)[i + 2], (*this)[i + 3], (*this)[i + 4], (*this)[i + 5]);
	}

	size_t avg_len() const
	{
		return this->letters() / this->get_length();
	}

	virtual ~Sequence_set()
	{ }

};
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

#ifndef REDUCTION_H_
#define REDUCTION_H_

#include <vector>
#include <string>
#include <string.h>
#include "value.h"
#include "../util/util.h"
#include "sequence.h"

using std::string;
using std::vector;

struct Reduction
{

	Reduction(const char *definition_string)
	{
		memset(map_, 0, sizeof(map_));
		memset(map8_, 0, sizeof(map8_));
		map_[(long)value_traits.mask_char] = value_traits.mask_char;
		const vector<string> tokens(tokenize(definition_string, " "));
		size_ = (unsigned)tokens.size();
		bit_size_ = (uint64_t)ceil(log(size_) / log(2));
		for (unsigned i = 0; i<size_; ++i)
			for (unsigned j = 0; j<tokens[i].length(); ++j) {
				const char ch = tokens[i][j];
				map_[(long)value_traits.from_char(ch)] = i;
				map8_[(long)value_traits.from_char(ch)] = i;
			}
		map8_[(long)value_traits.mask_char] = (char)size_;
	}

	unsigned size() const
	{
		return size_;
	}

	uint64_t bit_size() const
	{
		return bit_size_;
	}

	unsigned operator()(Letter a) const
	{
		return map_[(long)a];
	}

	const char* map8() const
	{
		return map8_;
	}

	inline friend std::ostream& operator<<(std::ostream &os, const Reduction &r)
	{
		for (unsigned i = 0; i < r.size_; ++i) {
			os << '[';
			for (unsigned j = 0; j < 20; ++j)
				if (r. map_[j] == i)
					os << value_traits.alphabet[j];
			os << ']';
		}
		return os;
	}

	static void reduce_seq(const sequence &seq, vector<char> &dst)
	{
		dst.clear();
		dst.resize(seq.length());
		for (unsigned i = 0; i < seq.length(); ++i)
			dst[i] = reduction(seq[i]);
	}

	static Reduction reduction;

private:

	unsigned map_[256];
#ifdef _MSC_VER
	__declspec(align(16)) char map8_[256];
#else
	char map8_[256] __attribute__((aligned(16)));
#endif
	unsigned size_;
	uint64_t bit_size_;

};

#ifdef EXTRA
#include "../../../extra/reduction.h"
#endif

#endif /* REDUCTION_H_ */
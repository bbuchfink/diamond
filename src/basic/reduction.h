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

};

#ifdef EXTRA
#include "../../../extra/reduction.h"
#endif

#endif /* REDUCTION_H_ */

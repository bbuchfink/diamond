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

using std::string;
using std::vector;

#include "value.h"

template<typename _val>
struct Reduction
{

    Reduction(const char *definition_string)
	{
		memset(map_, 0, sizeof(map_));
		memset(map8_, 0, sizeof(map8_));
		const vector<string> tokens (tokenize(definition_string, " "));
		size_ = tokens.size();
        for(unsigned i=0;i<size_;++i)
        	for(unsigned j=0;j<tokens[i].length();++j) {
        		const char ch = tokens[i][j];
        		map_[(long) Value_traits<_val>::from_char(ch)] =  i;
                map8_[(long) Value_traits<_val>::from_char(ch)] =  i;
            }
	}

	unsigned size() const
	{ return size_; }

	unsigned operator()(_val a) const
	{ return map_[(long)a]; }

	const char* map8() const
	{ return map8_; }

	static const Reduction reduction;

private:

	unsigned map_[256];
	char map8_[256];
	unsigned size_;

};

template<> const Reduction<Amino_acid> Reduction<Amino_acid>::reduction ("KREDQN C G H M F Y ILV W P STA");

#ifdef EXTRA
#include "../../extra/reduction.h"
#endif

#endif /* REDUCTION_H_ */

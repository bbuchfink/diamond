/****
Copyright (c) 2016, Benjamin Buchfink
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/

#ifndef INDEX_H_
#define INDEX_H_

#include <vector>
#include "sequence_set.h"
#include "../util/hash_table.h"

using std::vector;

struct Seed_double_index
{

	Seed_double_index()
	{}
	Seed_double_index(size_t psize);
	Seed_double_index(const Array<unsigned, Hashed_seed::p> &psize);
	
private:

	static const double hash_table_factor;
	
	struct Entry
	{
		uint32_t q;
		operator unsigned() const
		{
			return q;
		}
	};

	PHash_table<Entry> tables[Hashed_seed::p];

	friend struct Count_query_callback;
	friend struct Access;

};

extern Seed_double_index index[Const::max_shapes];

vector<Array<unsigned, Hashed_seed::p> > count_exact(const Sequence_set &seqs);
vector<size_t> count_approximate(const Sequence_set &seqs);
void build_query_index();

#endif
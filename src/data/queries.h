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

#ifndef QUERIES_H_
#define QUERIES_H_

#include "../basic/translate.h"
#include "../util/complexity_filter.h"
#include "sorted_list.h"
#include "../basic/statistics.h"
#include "sequence_set.h"

extern Partitioned_histogram query_hst;
extern unsigned current_query_chunk;

struct query_source_seqs
{
	static const Sequence_set& get()
	{ return *data_; }
	static Sequence_set *data_;
};

struct query_seqs
{
	static const Sequence_set& get()
	{ return *data_; }
	static Sequence_set *data_;
};

struct query_ids
{
	static const String_set<0>& get()
	{ return *data_; }
	static String_set<0> *data_;
};

extern vector<bool> query_aligned;

void write_unaligned(Output_stream *file);

inline unsigned get_source_query_len(unsigned query_id)
{
	return align_mode.query_translated ? (unsigned)query_seqs::get().reverse_translated_len(query_id*align_mode.query_contexts) : (unsigned)query_seqs::get().length(query_id);
}

#endif /* QUERIES_H_ */

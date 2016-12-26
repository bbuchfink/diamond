/****
Copyright (c) 2014-2016, University of Tuebingen, Benjamin Buchfink
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

#ifndef ALIGN_QUERIES_H_
#define ALIGN_QUERIES_H_

#include <deque>
#include "../util/merge_sort.h"
#include "../search/trace_pt_buffer.h"
#include "../util/map.h"
#include "../util/task_queue.h"
#include "align.h"
#include "query_mapper.h"

using std::vector;

struct Query_data
{
	Query_data():
		mapper(0)
	{}
	Query_data(Query_mapper *mapper):
		mapper(mapper),
		state(init)
	{}
	Query_mapper *mapper;
	Text_buffer buf;
	enum { init, free, closing, finished };
	unsigned state;
};

struct Query_queue
{
	void init(Trace_pt_list::iterator begin, Trace_pt_list::iterator end);
	void flush(Output_stream *out, Statistics &stat);
	Query_data* get();
	void pop_busy();

	tthread::mutex lock;
	std::deque<Query_data*> queue;
	std::queue<Query_data*> out_queue;
	Trace_pt_list::iterator trace_pt_pos, trace_pt_end;
	bool writing;
	unsigned n, last_query;
};

extern Query_queue query_queue;

void align_worker(Output_stream *out);
void align_queries(const Trace_pt_buffer &trace_pts, Output_stream* output_file);

#endif /* ALIGN_QUERIES_H_ */

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
void align_queries(Trace_pt_buffer &trace_pts, Output_stream* output_file);

#endif /* ALIGN_QUERIES_H_ */

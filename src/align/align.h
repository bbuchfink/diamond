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

#ifndef ALIGN_H_
#define ALIGN_H_

#include <memory>
#include <vector>
#include <map>
#include "../search/trace_pt_buffer.h"
#include "../util/task_queue.h"
#include "../basic/statistics.h"
#include "align_struct.h"

using std::vector;
using std::auto_ptr;

struct Output_writer
{
	Output_writer(Output_stream* f) :
		f_(f)
	{ }
	void operator()(Text_buffer &buf)
	{
		f_->write(buf.get_begin(), buf.size());
		buf.clear();
	}
private:
	Output_stream* const f_;
};

template<typename _buffer>
struct Ring_buffer_sink
{
	Ring_buffer_sink(Output_stream *output_file):
		writer(output_file),
		queue(config.threads_ * 4, writer)
	{}
	bool get(size_t &i, _buffer *& buffer, Trace_pt_list::Query_range &query_range)
	{
		return queue.get(i, buffer, query_range);
	}
	void push(size_t i)
	{
		queue.push(i);
	}
private:
	Output_writer writer;
	Task_queue<_buffer, Output_writer> queue;
};

void align_queries(Trace_pt_buffer &trace_pts, Output_stream* output_file);

#endif
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

#include <vector>
#include <map>
#include "../search/trace_pt_buffer.h"
#include "../util/task_queue.h"
#include "../basic/statistics.h"
#include "align_struct.h"

using std::vector;

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

struct Simple_query_queue
{
	enum { end = 0xffffffffffffffffllu };
	Simple_query_queue(size_t qbegin, size_t qend, vector<hit>::iterator begin, vector<hit>::iterator end):
		next_(qbegin),
		qend_(qend),
		it_(begin),
		end_(end)
	{}
	size_t get(vector<hit>::iterator &begin, vector<hit>::iterator &end);
	size_t next() const
	{
		return next_;
	}
	size_t qend() const
	{
		return qend_;
	}
	static Simple_query_queue& get()
	{
		return *instance;
	}
	static auto_ptr<Simple_query_queue> instance;
private:
	tthread::mutex mtx_;
	size_t next_;
	const size_t qend_;
	vector<hit>::iterator it_, end_;
};

struct Output_sink
{
	Output_sink(size_t begin, Output_stream *f):
		f_(f),
		next_(begin),
		size_(0),
		max_size_(0)
	{}
	void push(size_t n, Text_buffer *buf);
	size_t size() const
	{
		return size_;
	}
	size_t max_size() const
	{
		return max_size_;
	}
	static Output_sink& get()
	{
		return *instance;
	}
	size_t next() const
	{
		return next_;
	}
	static auto_ptr<Output_sink> instance;
private:
	void flush(Text_buffer *buf);
	tthread::mutex mtx_;
	Output_stream* const f_;
	std::map<size_t, Text_buffer*> backlog_;
	size_t next_, size_, max_size_;
};

#endif
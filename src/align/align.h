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

#ifndef ALIGN_H_
#define ALIGN_H_

#include <vector>
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

#endif
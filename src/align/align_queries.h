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

#ifndef ALIGN_QUERIES_H_
#define ALIGN_QUERIES_H_

#include "../util/merge_sort.h"
#include "../search/trace_pt_buffer.h"
#include "../util/map.h"
#include "align_read.h"
#include "../util/task_queue.h"
#include "align.h"

using std::vector;

template<unsigned _d>
void align_queries(Trace_pt_list::iterator begin,
		Trace_pt_list::iterator end,
		Output_buffer &buffer,
		Statistics &st)
{
	typedef Map<typename vector<hit>::iterator,typename hit::template Query_id<_d> > Map_t;
	Map_t hits (begin, end);
	typename Map_t::Iterator i = hits.begin();
	while(i.valid()) {
		align_read(buffer, st, i.begin(), i.end());
		++i;
	}
}

#define Output_sink Ring_buffer_sink

template<typename _buffer>
struct Align_context
{
	Align_context(Trace_pt_list &trace_pts, Output_stream* output_file):
		trace_pts (trace_pts),
		output_file (output_file),
		sink (output_file)
	{ }
	void operator()(unsigned thread_id)
	{
		Statistics st;
		size_t i=0;
		Trace_pt_list::Query_range query_range (trace_pts.get_range());
		_buffer *buffer = 0;
		while(sink.get(i, buffer, query_range)) {
			try {
				switch(query_contexts()) {
				case 6:
					align_queries<6>(query_range.begin, query_range.end, *buffer, st);
					break;
				case 2:
					align_queries<2>(query_range.begin, query_range.end, *buffer, st);
					break;
				default:
					align_queries<1>(query_range.begin, query_range.end, *buffer, st);
				}
				sink.push(i);
			}
			catch (std::bad_alloc&) {
				std::cout << "Out of memory error." << std::endl;
				std::terminate();
			}
			catch (std::exception &e) {
				std::cout << e.what() << std::endl;
				std::terminate();
				//queue.wake_all();
			}
		}
		statistics += st;
	}
	Trace_pt_list &trace_pts;
	Output_stream* output_file;
	Output_sink<_buffer> sink;
};

void align_queries(const Trace_pt_buffer &trace_pts, Output_stream* output_file)
{
	Trace_pt_list v;
	for(unsigned bin=0;bin<trace_pts.bins();++bin) {
		log_stream << "Processing query bin " << bin+1 << '/' << trace_pts.bins() << '\n';
		task_timer timer ("Loading trace points", 3);
		statistics.max(Statistics::TEMP_SPACE, trace_pts.load(v, bin));
		timer.go("Sorting trace points");
		merge_sort(v.begin(), v.end(), config.threads_);
		v.init();
		timer.go("Computing alignments");
		if(ref_header.n_blocks > 1) {
			Align_context<Temp_output_buffer> context (v, output_file);
			launch_thread_pool(context, config.threads_);
		} else {
			Align_context<Output_buffer> context (v, output_file);
			launch_thread_pool(context, config.threads_);
		}
	}
}

#endif /* ALIGN_QUERIES_H_ */

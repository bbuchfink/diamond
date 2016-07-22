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
#include "align_read.h"
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
	void init(Trace_pt_list::iterator begin, Trace_pt_list::iterator end)
	{
		trace_pt_pos = begin;
		trace_pt_end = end;
		assert(queue.empty());
		assert(out_queue.empty());
		n = 0;
	}
	void flush(Output_stream *out, Statistics &stat)
	{
		writing = true;
		std::queue<Query_data*> q;
		while(true) {
			while (!out_queue.empty() && out_queue.front()->state == Query_data::finished) {
				q.push(out_queue.front());
				out_queue.pop();
			}
			if (q.empty()) {
				writing = false;
				lock.unlock();
				return;
			}
			lock.unlock();

			unsigned k = 0;
			while (!q.empty()) {
				out->write(q.front()->buf.get_begin(), q.front()->buf.size());
				delete q.front();
				q.pop();
				++k;
			}

			lock.lock();
			n -= k;
			if (n > 100)
				cout << "qlen=" << out_queue.size() << " finished=" << n << endl;
		}
	}
	Query_data* get()
	{
		for (std::deque<Query_data*>::iterator i = queue.begin(); i != queue.end(); ++i)
			if ((*i)->state == Query_data::free)
				return *i;
		return 0;
	}
	void pop_busy()
	{
		while (!queue.empty() && (queue.front()->state == Query_data::closing || queue.front()->state == Query_data::finished)) {
			//cout << "pop " << queue.front()->query_id << endl;
			out_queue.push(queue.front());
			queue.pop_front();
		}
	}

	tthread::mutex lock;
	std::deque<Query_data*> queue;
	std::queue<Query_data*> out_queue;
	Trace_pt_list::iterator trace_pt_pos, trace_pt_end;
	bool writing;
	unsigned n;
};

extern Query_queue query_queue;

inline void align_worker(Output_stream *out)
{
	Statistics stat;
	Query_data *data = 0;
	unsigned n_targets = 0;
	
	while (true) {
		
		query_queue.lock.lock();

		if (data) {
			data->mapper->targets_finished += n_targets;			
			if (data->mapper->finished()) {
				query_queue.lock.unlock();
				data->mapper->generate_output(data->buf, stat);
				delete data->mapper;
				query_queue.lock.lock();
				data->state = Query_data::finished;
				++query_queue.n;
				if (!query_queue.writing && !query_queue.out_queue.empty() && data == query_queue.out_queue.front()) {
					query_queue.flush(out, stat);
					data = 0;
					continue;
				}
			}
		}
		
		if (!(data = query_queue.get())) {
			if (query_queue.trace_pt_pos >= query_queue.trace_pt_end) {
				//if (query_queue.queue.empty()) {
					query_queue.lock.unlock();
					//cout << "finished" << endl;
					break;
				//}			
				query_queue.lock.unlock();				
			}
			else {
				query_queue.queue.push_back(new Query_data(new Query_mapper()));
				data = query_queue.queue.back();
				//cout << "init " << mapper->query_id << " qlen=" << query_queue.queue.size() << endl;
				query_queue.lock.unlock();
				data->mapper->init();
				data->state = Query_data::free;
				data = 0;
			}
		}
		else {
			size_t target = data->mapper->next_target;
			n_targets = std::min(config.target_fetch_size, (unsigned)data->mapper->n_targets() - data->mapper->next_target);
			data->mapper->next_target += n_targets;
			if (target + n_targets == data->mapper->n_targets()) {
				data->state = Query_data::closing;
				query_queue.pop_busy();
			}
			query_queue.lock.unlock();

			for (unsigned i = 0; i < n_targets; ++i)
				data->mapper->align_target(target + i, stat);
		}

	}
	
	statistics += stat;
}

template<unsigned _d>
inline void align_queries(Trace_pt_list::iterator begin,
		Trace_pt_list::iterator end,
		Text_buffer &buffer,
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
		Text_buffer *buffer = 0;
		while(sink.get(i, buffer, query_range)) {
			try {
				switch(align_mode.query_contexts) {
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
			}
		}
		statistics += st;
	}
	Trace_pt_list &trace_pts;
	Output_stream* output_file;
	Output_sink<Text_buffer> sink;
};

inline void align_queries(const Trace_pt_buffer &trace_pts, Output_stream* output_file)
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
		if (config.target_parallel) {
			query_queue.init(v.begin(), v.end());
			Thread_pool threads;
			for (unsigned i = 0; i < config.threads_; ++i)
				threads.push_back(launch_thread(align_worker, output_file));
			threads.join_all();
		}
		else {
			Align_context context(v, output_file);
			launch_thread_pool(context, config.threads_);
		}
	}
}

#endif /* ALIGN_QUERIES_H_ */

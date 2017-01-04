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

#include "align_queries.h"
#include "../data/reference.h"
#include "../output/output_format.h"

Query_queue query_queue;

void Query_queue::init(Trace_pt_list::iterator begin, Trace_pt_list::iterator end)
{
	trace_pt_pos = begin;
	trace_pt_end = end;
	assert(queue.empty());
	assert(out_queue.empty());
	n = 0;
}
	
void Query_queue::flush(Output_stream *out, Statistics &stat)
{
	writing = true;
	std::queue<Query_data*> q;
	while (true) {
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
		/*if (n > 100)
		cout << "qlen=" << out_queue.size() << " finished=" << n << endl;*/
	}
}

Query_data* Query_queue::get()
{
	for (std::deque<Query_data*>::iterator i = queue.begin(); i != queue.end(); ++i)
		if ((*i)->state == Query_data::free)
			return *i;
	return 0;
}

void Query_queue::pop_busy()
{
	while (!queue.empty() && (queue.front()->state == Query_data::closing || queue.front()->state == Query_data::finished)) {
		out_queue.push(queue.front());
		queue.pop_front();
	}
}

void align_worker(Output_stream *out)
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
				const bool aligned = data->mapper->generate_output(data->buf, stat);
				const unsigned query_id = data->mapper->query_id;
				delete data->mapper;
				query_queue.lock.lock();
				data->state = Query_data::finished;
				if (aligned && !config.unaligned.empty())
					query_aligned[query_id] = true;
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
				query_queue.lock.unlock();
				break;
			}
			else {
				query_queue.queue.push_back(new Query_data(new Query_mapper()));
				data = query_queue.queue.back();
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

void align_queries(const Trace_pt_buffer &trace_pts, Output_stream* output_file)
{
	query_queue.last_query = (unsigned)-1;
	for (unsigned bin = 0; bin < trace_pts.bins(); ++bin) {
		log_stream << "Processing query bin " << bin + 1 << '/' << trace_pts.bins() << '\n';
		task_timer timer("Loading trace points", 3);
		Trace_pt_list *v = new Trace_pt_list;
		statistics.max(Statistics::TEMP_SPACE, trace_pts.load(*v, bin));
		timer.go("Sorting trace points");
		merge_sort(v->begin(), v->end(), config.threads_);
		v->init();
		timer.go("Computing alignments");
		query_queue.init(v->begin(), v->end());
		Thread_pool threads;
		for (unsigned i = 0; i < config.threads_; ++i)
			threads.push_back(launch_thread(align_worker, output_file));
		threads.join_all();
		timer.go("Deallocating buffers");
		delete v;
	}
	if (!blocked_processing && *output_format != Output_format::daa && config.report_unaligned != 0) {
		Text_buffer buf;
		for (unsigned i = query_queue.last_query + 1; i < query_ids::get().get_length(); ++i) {
			output_format->print_query_intro(i, query_ids::get()[i].c_str(), get_source_query_len(i), buf, true);
			output_format->print_query_epilog(buf, true);
		}
		output_file->write(buf.get_begin(), buf.size());
	}
}
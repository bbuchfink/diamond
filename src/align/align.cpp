/****
Copyright (c) 2016-2017, Benjamin Buchfink
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

#include "../basic/value.h"
#include "align.h"
#include "align_queries.h"
#include "../data/reference.h"
#include "../output/output_format.h"

using std::map;

auto_ptr<Simple_query_queue> Simple_query_queue::instance;
auto_ptr<Output_sink> Output_sink::instance;

size_t Simple_query_queue::get(vector<hit>::iterator &begin, vector<hit>::iterator &end)
{
	mtx_.lock();
	const unsigned query = (unsigned)(next_++),
		c = align_mode.query_contexts;
	if (query >= qend_) {
		mtx_.unlock();
		return Simple_query_queue::end;
	}
	begin = it_;
	while (it_ < end_ && it_->query_ / c == query)
		++it_;
	end = it_;
	mtx_.unlock();	
	return query;
}

void Output_sink::push(size_t n, Text_buffer *buf)
{
	mtx_.lock();
	//cout << "n=" << n << " next=" << next_ << endl;
	if (n != next_) {
		backlog_[n] = buf;
		size_ += buf ? buf->alloc_size() : 0;
		max_size_ = std::max(max_size_, size_);
		mtx_.unlock();
	}
	else
		flush(buf);
}

void Output_sink::flush(Text_buffer *buf)
{
	size_t n = next_ + 1;
	vector<Text_buffer*> out;
	out.push_back(buf);
	map<size_t, Text_buffer*>::iterator i;
	do {
		while ((i = backlog_.begin()) != backlog_.end() && i->first == n) {
			out.push_back(i->second);
			backlog_.erase(i);
			++n;
		}
		mtx_.unlock();
		size_t size = 0;
		for (vector<Text_buffer*>::iterator j = out.begin(); j < out.end(); ++j) {
			if (*j) {
				f_->write((*j)->get_begin(), (*j)->size());
				if(*j != buf)
					size += (*j)->alloc_size();
				delete *j;
			}
		}
		out.clear();
		mtx_.lock();
		size_ -= size;
	} while ((i = backlog_.begin()) != backlog_.end() && i->first == n);
	next_ = n;
	mtx_.unlock();
}

void align_worker(size_t thread_id)
{
	vector<hit>::iterator begin, end;
	size_t query;
	Statistics stat;
	while ((query = Simple_query_queue::get().get(begin, end)) != Simple_query_queue::end) {
		if (end == begin) {
			Output_sink::get().push(query, 0);
			continue;
		}
		Query_mapper mapper(query, begin, end);
		mapper.init();
		if (config.ext == Config::swipe)
			mapper.align_targets(stat);
		else {
			stat.inc(Statistics::TARGET_HITS0, mapper.n_targets());
			for (size_t i = 0; i < mapper.n_targets(); ++i)
				mapper.ungapped_stage(i);
			mapper.rank_targets(config.rank_ratio);
			stat.inc(Statistics::TARGET_HITS1, mapper.n_targets());
			for (size_t i = 0; i < mapper.n_targets(); ++i)
				mapper.greedy_stage(i);
			mapper.rank_targets(config.rank_ratio2);
			stat.inc(Statistics::TARGET_HITS2, mapper.n_targets());
			for (size_t i = 0; i < mapper.n_targets(); ++i)
				mapper.align_target(i, stat);
		}
		Text_buffer *buf = new Text_buffer;
		const bool aligned = mapper.generate_output(*buf, stat);
		if (aligned && !config.unaligned.empty())
			query_aligned[query] = true;
		Output_sink::get().push(query, buf);
	}
	statistics += stat;
}

void heartbeat_worker()
{
	static const int interval = 100;
	int n = 0;
	while (Output_sink::get().next() < Simple_query_queue::get().qend()) {
		if (n == interval) {
			const string title(query_ids::get()[Output_sink::get().next()].c_str());
			verbose_stream << "Queries=" << Simple_query_queue::get().next() << " size=" << megabytes(Output_sink::get().size()) << " max_size=" << megabytes(Output_sink::get().max_size())
				<< " next=" << title.substr(0, title.find(' ')) << endl;
			n = 0;
		}
		else
			++n;
		tthread::this_thread::sleep_for(tthread::chrono::milliseconds(10));
	}
}

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

void align_queries(Trace_pt_buffer &trace_pts, Output_stream* output_file)
{
	query_queue.last_query = (unsigned)-1;
	const size_t max_size = (size_t)std::min(config.chunk_size*1e9 * 9 * 2 / config.lowmem, 2e9);
	pair<size_t, size_t> query_range;
	while(true) {
		task_timer timer("Loading trace points", 3);
		Trace_pt_list *v = new Trace_pt_list;
		statistics.max(Statistics::TEMP_SPACE, trace_pts.load(*v, max_size, query_range));
		if (query_range.second - query_range.first == 0) {
			delete v;
			break;
		}
		timer.go("Sorting trace points");
		merge_sort(v->begin(), v->end(), config.threads_);
		v->init();
		timer.go("Computing alignments");
		if (config.load_balancing == Config::target_parallel) {
			query_queue.init(v->begin(), v->end());
			Thread_pool threads;
			for (unsigned i = 0; i < config.threads_; ++i)
				threads.push_back(launch_thread(static_cast<void (*)(Output_stream*)>(&align_worker), output_file));
			threads.join_all();
		}
		else {
			Simple_query_queue::instance = auto_ptr<Simple_query_queue>(new Simple_query_queue(query_range.first, query_range.second, v->begin(), v->end()));
			Output_sink::instance = auto_ptr<Output_sink>(new Output_sink(query_range.first, output_file));
			Thread_pool threads;
			if (config.verbosity >= 3)
				threads.push_back(launch_thread(heartbeat_worker));
			for (size_t i = 0; i < config.threads_; ++i)
				threads.push_back(launch_thread(static_cast<void(*)(size_t)>(&align_worker), i));
			threads.join_all();
		}
		timer.go("Deallocating buffers");
		delete v;
	}
	if (!blocked_processing && *output_format != Output_format::daa && config.report_unaligned != 0 && config.load_balancing == Config::target_parallel) {
		Text_buffer buf;
		for (unsigned i = query_queue.last_query + 1; i < query_ids::get().get_length(); ++i) {
			output_format->print_query_intro(i, query_ids::get()[i].c_str(), get_source_query_len(i), buf, true);
			output_format->print_query_epilog(buf, true);
		}
		output_file->write(buf.get_begin(), buf.size());
	}
}
/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#include <memory>
#include "../basic/value.h"
#include "align.h"
#include "../data/reference.h"
#include "../output/output_format.h"
#include "../util/queue.h"
#include "../output/output.h"
#include "query_mapper.h"
#include "../util/merge_sort.h"

using namespace std;

DpStat dp_stat;

struct Align_fetcher
{
	static void init(size_t qbegin, size_t qend, vector<hit>::iterator begin, vector<hit>::iterator end)
	{
		it_ = begin;
		end_ = end;
		queue_ = unique_ptr<Queue>(new Queue(qbegin, qend));
	}
	void operator()(size_t query)
	{
		const unsigned q = (unsigned)query,
			c = align_mode.query_contexts;
		begin = it_;
		while (it_ < end_ && it_->query_ / c == q)
			++it_;
		end = it_;
		this->query = query;
	}
	bool get()
	{
		return queue_->get(*this) != Queue::end;
	}
	size_t query;
	vector<hit>::iterator begin, end;
private:	
	static vector<hit>::iterator it_, end_;
	static unique_ptr<Queue> queue_;
};

unique_ptr<Queue> Align_fetcher::queue_;
vector<hit>::iterator Align_fetcher::it_;
vector<hit>::iterator Align_fetcher::end_;

void align_worker(size_t thread_id, const Parameters *params, const Metadata *metadata)
{
	Align_fetcher hits;
	Statistics stat;
	DpStat dp_stat;
	while (hits.get()) {
		if (hits.end == hits.begin) {
			TextBuffer *buf = 0;
			if (!blocked_processing && *output_format != Output_format::daa && config.report_unaligned != 0) {
				buf = new TextBuffer;
				const char *query_title = query_ids::get()[hits.query].c_str();
				output_format->print_query_intro(hits.query, query_title, get_source_query_len((unsigned)hits.query), *buf, true);
				output_format->print_query_epilog(*buf, query_title, true, *params);
			}
			OutputSink::get().push(hits.query, buf);
			continue;
		}

		QueryMapper *mapper;
		if (config.ext == Config::swipe)
			mapper = new ExtensionPipeline::Swipe::Pipeline(*params, hits.query, hits.begin, hits.end);
		else if (config.frame_shift != 0)
			mapper = new ExtensionPipeline::BandedSwipe::Pipeline(*params, hits.query, hits.begin, hits.end, dp_stat);
		else
			mapper = new ExtensionPipeline::Greedy::Pipeline(*params, hits.query, hits.begin, hits.end);
		task_timer timer("Initializing mapper", config.load_balancing == Config::target_parallel ? 3 : UINT_MAX);
		mapper->init();
		timer.finish();
		mapper->run(stat);

		timer.go("Generating output");
		TextBuffer *buf = 0;
		if (*output_format != Output_format::null) {
			buf = new TextBuffer;
			const bool aligned = mapper->generate_output(*buf, stat, *metadata);
			if (aligned && (!config.unaligned.empty() || !config.aligned_file.empty()))
				query_aligned[hits.query] = true;
		}
		delete mapper;
		OutputSink::get().push(hits.query, buf);
	}
	statistics += stat;
	::dp_stat += dp_stat;
}

void align_queries(Trace_pt_buffer &trace_pts, Consumer* output_file, const Parameters &params, const Metadata &metadata)
{
	const size_t max_size = (size_t)std::min(config.chunk_size*1e9 * 9 * 2 / config.lowmem, 2e9);
	pair<size_t, size_t> query_range;
	while (true) {
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
		Align_fetcher::init(query_range.first, query_range.second, v->begin(), v->end());
		OutputSink::instance = unique_ptr<OutputSink>(new OutputSink(query_range.first, output_file));
		Thread_pool threads;
		if (config.verbosity >= 3 && config.load_balancing == Config::query_parallel)
			threads.push_back(launch_thread(heartbeat_worker, query_range.second));
		size_t n_threads = config.load_balancing == Config::query_parallel ? (config.threads_align == 0 ? config.threads_ : config.threads_align) : 1;
		for (size_t i = 0; i < n_threads; ++i)
			//threads.push_back(launch_thread(static_cast<void(*)(size_t)>(&align_worker), i));
			threads.push_back(launch_thread(align_worker, i, &params, &metadata));
		threads.join_all();
		timer.finish();

		double t = timer.get();
		log_stream << "Gross cells = " << dp_stat.gross_cells << endl;
		log_stream << "Gross GCUPS = " << (double)dp_stat.gross_cells / 1e9 / t << endl;
		log_stream << "Gross GCUPS/thread = " << (double)dp_stat.gross_cells / n_threads / 1e9 / t << endl;
		log_stream << "Net cells = " << dp_stat.net_cells << endl;
		log_stream << "Net GCUPS = " << (double)dp_stat.net_cells / 1e9 / t << endl;
		log_stream << "Net GCUPS/thread = " << (double)dp_stat.net_cells / n_threads / 1e9 / t << endl;

		timer.go("Deallocating buffers");
		delete v;
	}
}
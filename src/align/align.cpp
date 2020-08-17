/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen

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
#include "legacy/query_mapper.h"
#include "../util/merge_sort.h"
#include "extend.h"
#include "../util/algo/radix_sort.h"

using std::get;
using std::tuple;
using std::unique_ptr;

DpStat dp_stat;

struct Align_fetcher
{
	static void init(size_t qbegin, size_t qend, hit* begin, hit* end)
	{
		it_ = begin;
		end_ = end;
		queue_ = unique_ptr<Queue>(new Queue(qbegin, qend));
	}
	bool operator()(size_t query)
	{
		const unsigned q = (unsigned)query,
			c = align_mode.query_contexts;
		begin = it_;
		while (it_ < end_ && it_->query_ / c == q)
			++it_;
		end = it_;
		this->query = query;
		target_parallel = (end - begin > config.query_parallel_limit) && (config.frame_shift == 0 || (config.toppercent < 100 && config.query_range_culling));
		return target_parallel;
	}
	bool get()
	{
		return queue_->get(*this) != Queue::end;
	}
	void release() {
		if (target_parallel)
			queue_->release();
	}
	size_t query;
	hit* begin, *end;
	bool target_parallel;
private:	
	static hit* it_, *end_;
	static unique_ptr<Queue> queue_;
};

unique_ptr<Queue> Align_fetcher::queue_;
hit* Align_fetcher::it_;
hit* Align_fetcher::end_;

TextBuffer* legacy_pipeline(Align_fetcher &hits, const sequence *subjects, size_t subject_count, const Metadata *metadata, const Parameters *params, Statistics &stat) {
	if ((hits.end == hits.begin) && subjects == nullptr) {
		TextBuffer *buf = nullptr;
		if (!blocked_processing && *output_format != Output_format::daa && config.report_unaligned != 0) {
			buf = new TextBuffer;
			const char *query_title = query_ids::get()[hits.query];
			output_format->print_query_intro(hits.query, query_title, get_source_query_len((unsigned)hits.query), *buf, true);
			output_format->print_query_epilog(*buf, query_title, true, *params);
		}
		return buf;
	}

	QueryMapper *mapper = new ExtensionPipeline::BandedSwipe::Pipeline(*params, hits.query, hits.begin, hits.end, dp_stat, *metadata, hits.target_parallel);

	task_timer timer("Initializing mapper", hits.target_parallel ? 3 : UINT_MAX);
	mapper->init();
	timer.finish();
	mapper->run(stat, subjects, subject_count);

	timer.go("Generating output");
	TextBuffer *buf = nullptr;
	if (*output_format != Output_format::null) {
		buf = new TextBuffer;
		const bool aligned = mapper->generate_output(*buf, stat);
		if (aligned && (!config.unaligned.empty() || !config.aligned_file.empty())) {
			query_aligned_mtx.lock();
			query_aligned[hits.query] = true;
			query_aligned_mtx.unlock();
		}
	}
	delete mapper;
	return buf;
}

void align_worker(size_t thread_id, const Parameters *params, const Metadata *metadata, const sequence *subjects, size_t subject_count)
{
	Align_fetcher hits;
	Statistics stat;
	DpStat dp_stat;
	while (hits.get()) {
		if(config.frame_shift != 0) {
			TextBuffer *buf = legacy_pipeline(hits, subjects, subject_count, metadata, params, stat);
			OutputSink::get().push(hits.query, buf);
			hits.release();
			continue;
		}
		task_timer timer;
		vector<Extension::Match> matches = Extension::extend(*params, hits.query, hits.begin, hits.end, *metadata, stat, hits.target_parallel ? Extension::TARGET_PARALLEL : 0);
		TextBuffer *buf = Extension::generate_output(matches, hits.query, stat, *metadata, *params);
		if (!matches.empty() && (!config.unaligned.empty() || !config.aligned_file.empty())) {
			query_aligned_mtx.lock();
			query_aligned[hits.query] = true;
			query_aligned_mtx.unlock();
		}
		OutputSink::get().push(hits.query, buf);
		if (hits.target_parallel)
			stat.inc(Statistics::TIME_TARGET_PARALLEL, timer.microseconds());
		hits.release();
	}
	statistics += stat;
	::dp_stat += dp_stat;
}

void align_queries(Trace_pt_buffer &trace_pts, Consumer* output_file, const Parameters &params, const Metadata &metadata)
{
	size_t max_size = std::min(size_t(config.chunk_size*1e9 * 10 * 2) / config.lowmem / 3, config.trace_pt_fetch_size);
	if (config.memory_limit != 0.0)
		max_size = std::max(max_size, size_t(config.memory_limit * 1e9));
	pair<size_t, size_t> query_range;
	vector<sequence> subjects;
	if (config.swipe_all) {
		subjects.reserve(ref_seqs::get().get_length());
		for (size_t i = 0; i < ref_seqs::get().get_length(); ++i)
			subjects.push_back(ref_seqs::get()[i]);
	}

	trace_pts.load(max_size);

	while (true) {
		task_timer timer("Loading trace points", 3);
		tuple<vector<hit>*, size_t, size_t> input = trace_pts.retrieve();
		if (get<0>(input) == nullptr)
			break;
		statistics.inc(Statistics::TIME_LOAD_SEED_HITS, timer.microseconds());
		vector<hit>* hit_buf = get<0>(input);
		query_range = { get<1>(input), get<2>(input) };
		trace_pts.load(max_size);

		timer.go("Sorting trace points");
		//if (config.beta)
			radix_sort<hit, hit::Query>(hit_buf->data(), hit_buf->data() + hit_buf->size(), (uint32_t)query_range.second * align_mode.query_contexts, config.threads_);
		//else
			//merge_sort(hit_buf->begin(), hit_buf->end(), config.threads_);
		statistics.inc(Statistics::TIME_SORT_SEED_HITS, timer.microseconds());

		timer.go("Computing alignments");
		Align_fetcher::init(query_range.first, query_range.second, hit_buf->data(), hit_buf->data() + hit_buf->size());
		OutputSink::instance = unique_ptr<OutputSink>(new OutputSink(query_range.first, output_file));
		vector<std::thread> threads;
		if (config.verbosity >= 3 && config.load_balancing == Config::query_parallel && !config.no_heartbeat)
			threads.emplace_back(heartbeat_worker, query_range.second);
		size_t n_threads = config.load_balancing == Config::query_parallel ? (config.threads_align == 0 ? config.threads_ : config.threads_align) : 1;
		for (size_t i = 0; i < n_threads; ++i)
			threads.emplace_back(align_worker, i, &params, &metadata, subjects.empty() ? nullptr : subjects.data(), subjects.size());
		for (auto &t : threads)
			t.join();
		statistics.inc(Statistics::TIME_EXT, timer.microseconds());
		
		timer.go("Deallocating buffers");
		delete hit_buf;
	}
	statistics.max(Statistics::SEARCH_TEMP_SPACE, trace_pts.total_disk_size());
}
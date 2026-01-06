/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

#ifdef _WIN32
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif
#include <memory>
#include "basic/value.h"
#include "align.h"
#include "output/output_format.h"
#include "output/output.h"
#include "legacy/pipeline.h"
#include "search/hit_buffer.h"
#include "util/parallel/thread_pool.h"
#include "extend.h"
#include "util/util.h"
#ifdef WITH_DNA
#include "../dna/extension.h"
#endif
#include "util/util.h"
#define _REENTRANT
#include "ips4o/ips4o.hpp"
#include "data/queries.h"
#include "data/sequence_file.h"
#include "search/hit_buffer.h"

using std::get;
using std::tuple;
using std::unique_ptr;
using std::thread;
using std::lock_guard;
using std::mutex;
using std::pair;
using std::vector;

DpStat dp_stat;

static vector<int64_t> make_partition(Search::Hit* begin, Search::Hit* end) {
	vector<int64_t> partition;
	partition.reserve(div_up(end - begin, (ptrdiff_t)config.min_task_trace_pts) + 1);
	Search::Hit* p = begin;
	partition.push_back(0);
	const BlockId c = align_mode.query_contexts;
	while (p < end) {
		Search::Hit* q = std::min(p + config.min_task_trace_pts, end - 1);
		const BlockId query = q->query_ / c;
		do {
			++q;
		} while (q < end && q->query_ / c == query);
		partition.push_back(q - begin);
		p = q;
	}
	return partition;
}

struct HitIterator {
	static bool single_query() {
		return config.swipe_all || align_mode.mode == AlignMode::blastn;
	}
	HitIterator(BlockId qbegin, BlockId qend, Search::Hit* begin, Search::Hit* end, vector<int64_t>::const_iterator partition, int64_t parts):
		partition(partition),
		parts(parts),
		data(begin),
		query_begin(qbegin),
		query_end(qend)
	{}
	struct Hits {
		BlockId query;
		Search::Hit* begin, * end;
	};
	vector<Hits> fetch(int64_t i) {
		vector<Hits> r;
		if (single_query()) {
			r.push_back(Hits{ (BlockId)i,nullptr,nullptr });
			return r;
		}
		assert(i >= 0 && i < parts);
		const BlockId c = align_mode.query_contexts;
		Search::Hit* begin = data + partition[i], * end = data + partition[i + 1];
		BlockId last_query = begin > data ? (begin - 1)->query_/c + 1 : query_begin;
		const int64_t query_count = (end - 1)->query_/c + 1 - last_query;
		r.reserve(query_count);
		for (; last_query < begin->query_/c; ++last_query)
			r.push_back(Hits{ last_query, nullptr,nullptr });
		auto it = merge_keys(begin, end, [c](const Search::Hit& h) { return h.query_/c; });
		while (it.good()) {
			for (; last_query < it.key(); ++last_query)
				r.push_back(Hits{ last_query, nullptr,nullptr });
			r.push_back(Hits{ (BlockId)it.key(), it.begin(), it.end() });
			++it;
			++last_query;
		}
		if (i == parts - 1) {
			r.reserve(r.size() + query_end - (r.back().query + 1));
			for (BlockId j = r.back().query + 1; j < query_end; ++j)
				r.push_back(Hits{ j, nullptr,nullptr });
		}
		return r;
	}
	const vector<int64_t>::const_iterator partition;
	const int64_t parts;
	Search::Hit* data;
	const BlockId query_begin, query_end;
};

static TextBuffer* legacy_pipeline(const HitIterator::Hits& hits, Search::Config& cfg, Statistics &stat) {
	if (hits.end == hits.begin) {
		TextBuffer *buf = nullptr;
		if (!cfg.blocked_processing && *cfg.output_format != OutputFormat::daa && cfg.output_format->report_unaligned()) {
			buf = new TextBuffer;
			Output::Info info{ cfg.query->seq_info(hits.query), true, cfg.db.get(), *buf, Util::Seq::AccessionParsing(), cfg.db->sequence_count(), cfg.db->letters() };
			cfg.output_format->print_query_intro(info);
			cfg.output_format->print_query_epilog(info);
		}
		return buf;
	}

	QueryMapper *mapper = new ExtensionPipeline::BandedSwipe::Pipeline(hits.query, hits.begin, hits.end, dp_stat, cfg);

	TaskTimer timer("Initializing mapper", UINT_MAX);
	mapper->init();
	timer.finish();
	mapper->run(stat, cfg);

	timer.go("Generating output");
	TextBuffer *buf = nullptr;
	if (*cfg.output_format != OutputFormat::null) {
		buf = new TextBuffer;
		const bool aligned = mapper->generate_output(*buf, stat, cfg);
		if (aligned && cfg.track_aligned_queries) {
			query_aligned_mtx.lock();
			if (!query_aligned[hits.query]) {
				query_aligned[hits.query] = true;
				++cfg.iteration_query_aligned;
			}
			query_aligned_mtx.unlock();
		}
	}
	delete mapper;
	return buf;
}

static void align_worker(HitIterator* hit_it, Search::Config* cfg, int64_t next)
{
	try {
		std::pmr::monotonic_buffer_resource pool;
		const vector<HitIterator::Hits> hits = hit_it->fetch(next);
		assert(!hits.empty());
		Statistics stat;
		DpStat dp_stat;
		const bool parallel = config.swipe_all && (cfg->target->seqs().size() >= cfg->query->seqs().size());

		for (auto h = hits.cbegin(); h < hits.cend(); ++h) {
			if (config.frame_shift != 0) {
				TextBuffer* buf = legacy_pipeline(*h, *cfg, stat);
				output_sink->push(h->query, buf);
				continue;
			}
			if (h->begin == nullptr && !HitIterator::single_query()) {
				output_sink->push(h->query, nullptr);
				continue;
			}

			vector<Extension::Match> matches =
#ifdef WITH_DNA
				align_mode.mode == AlignMode::blastn ? Dna::extend(*cfg, cfg->query->seqs()[h->query]) :
#endif
				Extension::extend(h->query, h->begin, h->end, *cfg, stat, parallel ? DP::Flags::PARALLEL : DP::Flags::NONE, pool);
			TextBuffer* buf = cfg->blocked_processing ? Extension::generate_intermediate_output(matches, h->query, *cfg) : Extension::generate_output(matches, h->query, stat, *cfg);
			if (!matches.empty() && cfg->track_aligned_queries) {
				std::lock_guard<std::mutex> lock(query_aligned_mtx);
				if (!query_aligned[h->query]) {
					query_aligned[h->query] = true;
					++cfg->iteration_query_aligned;
				}
			}
			if (!config.unaligned_targets.empty()) {
				lock_guard<mutex> lock(cfg->aligned_targets_mtx);
				for (const Extension::Match &m : matches) {
					const OId oid = cfg->target->block_id2oid(m.target_block_id);
					if (!cfg->aligned_targets[oid])
						cfg->aligned_targets[oid] = true;
				}
			}
			output_sink->push(h->query, buf);
		}

		statistics += stat;
		::dp_stat += dp_stat;
	}
	catch (std::exception& e) {
		exit_with_error(e);
	}
}

void align_queries(Consumer* output_file, Search::Config& cfg)
{
	const int64_t mem_limit = Util::String::interpret_number(config.memory_limit.get("16G"));

	pair<BlockId, BlockId> query_range;
	TaskTimer timer("Allocating memory", 3);

	if (!cfg.blocked_processing && !cfg.iterated())
		cfg.db->init_random_access(cfg.current_query_block, 0, false);

	int64_t res_size = cfg.query->mem_size() + cfg.target->mem_size(), last_size = 0;
	cfg.seed_hit_buf->alloc_buffer();
	//cfg.seed_hit_buf->load(std::min(mem_limit - res_size - cfg.seed_hit_buf->bin_size(1) * (int64_t)sizeof(Search::Hit), config.trace_pt_fetch_size));

	while (true) {
		timer.go("Loading trace points");
		if (!cfg.seed_hit_buf->load(std::min(mem_limit - res_size - cfg.seed_hit_buf->bin_size(1) * (int64_t)sizeof(Search::Hit), config.trace_pt_fetch_size)))
			break;
		tuple<Search::Hit*, int64_t, BlockId, BlockId> input = cfg.seed_hit_buf->retrieve();
		statistics.inc(Statistics::TIME_LOAD_SEED_HITS, timer.microseconds());
		timer.finish();
		Search::Hit* hit_buf = get<0>(input);
		const int64_t hit_count = get<1>(input);
		log_stream << "Processing " << hit_count << " trace points (" << Util::String::format(int64_t(hit_count * sizeof(Search::Hit))) << ")." << std::endl;
		res_size += hit_count * sizeof(Search::Hit);
		query_range = { get<2>(input), get<3>(input) };
		//cfg.seed_hit_buf->load(std::min(mem_limit - res_size, config.trace_pt_fetch_size));

		if (res_size + last_size > mem_limit)
			log_stream << "Warning: resident size (" << (res_size + last_size) << ") exceeds memory limit." << std::endl;

		timer.go("Sorting trace points");
#ifdef NDEBUG
		//sort::sort_parallel_blocked_inplace(hit_buf, hit_buf + hit_count, std::less<Search::Hit>(), config.threads_);
		ips4o::parallel::sort(hit_buf, hit_buf + hit_count, std::less<Search::Hit>(), config.threads_);
#else
		std::sort(hit_buf, hit_buf + hit_count);
#endif
		statistics.inc(Statistics::TIME_SORT_SEED_HITS, timer.microseconds());

		timer.go("Computing partition");
		const vector<int64_t> partition = make_partition(hit_buf, hit_buf + hit_count);

		timer.go("Computing alignments");
		HitIterator hit_it(query_range.first, query_range.second, hit_buf, hit_buf + hit_count, partition.begin(), (int64_t)partition.size() - 1);
        OutputWriter writer{output_file, cfg.blocked_processing ? '\0' : cfg.output_format->query_separator };
		output_sink.reset(new ReorderQueue<TextBuffer*, OutputWriter>(query_range.first, writer));
		unique_ptr<thread> heartbeat;
		if (config.verbosity >= 3 && config.load_balancing == Config::query_parallel && !config.swipe_all && config.heartbeat)
			heartbeat.reset(new thread(heartbeat_worker, query_range.second, &cfg));
		const int threads = config.load_balancing == Config::target_parallel || (config.swipe_all && (cfg.target->seqs().size() >= cfg.query->seqs().size())) ? 1
			: (config.threads_align == 0 ? config.threads_ : config.threads_align);
		auto task = [&hit_it, &cfg](ThreadPool& tp, int64_t i) { return align_worker(&hit_it, &cfg, i); };
		cfg.thread_pool.reset(config.swipe_all ? new ThreadPool(task, query_range.first, query_range.second) : new ThreadPool(task, 0, (int64_t)partition.size() - 1));
		cfg.thread_pool->run(threads);
		cfg.thread_pool->join();
		if (heartbeat)
			heartbeat->join();
		statistics.inc(Statistics::TIME_EXT, timer.microseconds());
		
		timer.go("Deallocating buffers");
		cfg.thread_pool.reset();
		output_sink.reset();
		last_size = hit_count * sizeof(Search::Hit);
		res_size -= last_size;
	}
	statistics.max(Statistics::SEARCH_TEMP_SPACE, cfg.seed_hit_buf->total_disk_size());

	timer.go("Freeing memory");
	cfg.seed_hit_buf->free_buffer();
	if (!cfg.blocked_processing && !cfg.iterated())
		cfg.db->end_random_access(false);
}
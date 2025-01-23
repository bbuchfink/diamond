/****
DIAMOND protein aligner
Copyright (C) 2013-2021 Max Planck Society for the Advancement of Science e.V.
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
#include "basic/value.h"
#include "align.h"
#include "output/output_format.h"
#include "output/output.h"
#include "legacy/pipeline.h"
#include "util/async_buffer.h"
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

using std::get;
using std::tuple;
using std::unique_ptr;
using std::thread;
using std::lock_guard;
using std::mutex;
using std::pair;
using std::vector;

DpStat dp_stat;

namespace Extension {

TextBuffer* pipeline_short(BlockId query, Search::Hit* begin, Search::Hit* end, Search::Config& cfg, Statistics& stats);

}

static vector<int64_t> make_partition(Search::Hit* begin, Search::Hit* end) {
	vector<int64_t> partition;
	partition.reserve(div_up(end - begin, config.min_task_trace_pts) + 1);
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
		if (!cfg.blocked_processing && *cfg.output_format != OutputFormat::daa && config.report_unaligned != 0) {
			buf = new TextBuffer;
			Output::Info info{ cfg.query->seq_info(hits.query), true, cfg.db.get(), *buf, {}, Util::Seq::AccessionParsing(), cfg.db->sequence_count(), cfg.db->letters() };
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
			if (config.pipeline_short) {
				TextBuffer* buf = Extension::pipeline_short(h->query, h->begin, h->end, *cfg, stat);
				output_sink->push(h->query, buf);
				continue;
			}
			if (h->begin == nullptr && !HitIterator::single_query()) {
				output_sink->push(h->query, nullptr);
				continue;
			}

			pair<vector<Extension::Match>, Extension::Stats> matches =
#ifdef WITH_DNA
				align_mode.mode == AlignMode::blastn ? Dna::extend(*cfg, cfg->query->seqs()[h->query]) :
#endif
				Extension::extend(h->query, h->begin, h->end, *cfg, stat, parallel ? DP::Flags::PARALLEL : DP::Flags::NONE);
			TextBuffer* buf = cfg->blocked_processing ? Extension::generate_intermediate_output(matches.first, h->query, *cfg) : Extension::generate_output(matches.first, matches.second, h->query, stat, *cfg);
			if (!matches.first.empty() && cfg->track_aligned_queries) {
				std::lock_guard<std::mutex> lock(query_aligned_mtx);
				if (!query_aligned[h->query]) {
					query_aligned[h->query] = true;
					++cfg->iteration_query_aligned;
				}
			}
			if (!config.unaligned_targets.empty()) {
				lock_guard<mutex> lock(cfg->aligned_targets_mtx);
				for (const Extension::Match &m : matches.first) {
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
	TaskTimer timer(nullptr, 3);

	if (!cfg.blocked_processing && !cfg.iterated())
		cfg.db->init_random_access(cfg.current_query_block, 0, false);

	int64_t res_size = cfg.query->mem_size() + cfg.target->mem_size(), last_size = 0;
	cfg.seed_hit_buf->load(std::min(mem_limit - res_size - cfg.seed_hit_buf->bin_size(1) * (int64_t)sizeof(Search::Hit), config.trace_pt_fetch_size));

	while (true) {
		timer.go("Loading trace points");
		tuple<vector<Search::Hit>*, BlockId, BlockId> input = cfg.seed_hit_buf->retrieve();
		if (get<0>(input) == nullptr)
			break;
		statistics.inc(Statistics::TIME_LOAD_SEED_HITS, timer.microseconds());
		vector<Search::Hit>* hit_buf = get<0>(input);
		res_size += hit_buf->size() * sizeof(Search::Hit);
		query_range = { get<1>(input), get<2>(input) };
		cfg.seed_hit_buf->load(std::min(mem_limit - res_size, config.trace_pt_fetch_size));

		if (res_size + last_size > mem_limit)
			log_stream << "Warning: resident size (" << (res_size + last_size) << ") exceeds memory limit." << std::endl;

		timer.go("Sorting trace points");
#ifdef NDEBUG
		ips4o::parallel::sort(hit_buf->begin(), hit_buf->end(), std::less<Search::Hit>(), config.threads_);
#else
		std::sort(hit_buf->begin(), hit_buf->end());
#endif
		statistics.inc(Statistics::TIME_SORT_SEED_HITS, timer.microseconds());

		timer.go("Computing partition");
		const vector<int64_t> partition = make_partition(hit_buf->data(), hit_buf->data() + hit_buf->size());

		timer.go("Computing alignments");
		HitIterator hit_it(query_range.first, query_range.second, hit_buf->data(), hit_buf->data() + hit_buf->size(), partition.begin(), (int64_t)partition.size() - 1);
        OutputWriter writer{output_file, cfg.output_format->query_separator};
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
		last_size = hit_buf->size() * sizeof(Search::Hit);
		res_size -= last_size;
		delete hit_buf;
	}
	statistics.max(Statistics::SEARCH_TEMP_SPACE, cfg.seed_hit_buf->total_disk_size());

	if (!cfg.blocked_processing && !cfg.iterated())
		cfg.db->end_random_access(false);
}
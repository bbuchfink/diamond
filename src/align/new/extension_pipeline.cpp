/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

#include <cstdint>
#include <string>
#include <algorithm>
#include "basic/match.h"
#include "basic/statistics.h"
#include "output/output_format.h"
#include "data/block/block.h"
#include "data/sequence_file.h"
#include "output/output.h"
#include "align/target.h"
#include "sw.h"
#include "thread_pool.h"
#include "extension_pipeline.h"

namespace ExtensionPipeline {

QueryState::QueryState(Run* run, BlockId query_id, Search::Hit* hits_end, Statistics& stats) :
	run(run),
	query_id(query_id),
	hits_end(hits_end),
	query(query_id, stats, run->cfg, *std::pmr::get_default_resource()),
	refs(1),
	hit_num(0),
	out(new TextBuffer)
{
}

Run::Run(Search::Hit* hits, uint64_t hit_count, const Search::Config& cfg):
	hits(hits),
	hit_count(hit_count),
	cfg(cfg),
	output_format(cfg.output_format->clone()),
	queries_pending(1),
	all_queries_submitted(false),
	pool(config.threads_, [this](Extension& e) { this->process_extension(e); })
{
}

static void extend(QueryState* state, BlockId target_id, const Search::Config& cfg, OutputFormat* format) {
	const BlockId query_id = state->query_id;
	const Sequence& query = cfg.query->seqs()[query_id];
	const Sequence& target = cfg.target->seqs()[target_id];
	Hsp hsp;
	if (!smith_waterman(query, target, hsp))
		return;

	const OId target_oid = cfg.target->block_id2oid(target_id);
	const bool all_seqids = flag_any(format->flags, Output::Flags::ALL_SEQIDS);
	const std::string target_title = cfg.target->has_ids()
		? cfg.target->ids()[target_id]
		: (flag_any(format->flags, Output::Flags::SSEQID) ? cfg.db->seqid(target_oid, all_seqids, true) : "");
	const Sequence& output_target = cfg.target->unmasked_seqs().empty() ? target : cfg.target->unmasked_seqs()[target_id];
	const double query_self_score = cfg.query->has_self_aln() ? cfg.query->self_aln_score(query_id) : 0.0;
	const double target_self_score = cfg.target->has_self_aln() ? cfg.target->self_aln_score(target_id) : 0.0;

	std::lock_guard<std::mutex> lock(state->mtx);
	const unsigned hit_num = state->hit_num++;
	Output::Info info{ cfg.query->seq_info(query_id), false, cfg.db.get(), *state->out, Util::Seq::AccessionParsing(), cfg.db->sequence_count(), cfg.db->letters() };
	format->print_match(HspContext(hsp,
		query_id,
		cfg.query->block_id2oid(query_id),
		TranslatedSequence(query),
		cfg.query->ids()[query_id],
		target_oid,
		static_cast<unsigned>(target.length()),
		target_title.c_str(),
		hit_num,
		0,
		output_target,
		0,
		query_self_score,
		target_self_score), info);
	statistics.inc(Statistics::MATCHES);
	statistics.inc(Statistics::PAIRWISE);
}

static void print_query_intro(BlockId query_id, const Search::Config& cfg, const OutputFormat& format, TextBuffer& out) {
	Output::Info info{ cfg.query->seq_info(query_id), false, cfg.db.get(), out, Util::Seq::AccessionParsing(), cfg.db->sequence_count(), cfg.db->letters() };
	format.print_query_intro(info);
}

static void print_query_epilog(BlockId query_id, const Search::Config& cfg, const OutputFormat& format, TextBuffer& out) {
	Output::Info info{ cfg.query->seq_info(query_id), false, cfg.db.get(), out, Util::Seq::AccessionParsing(), cfg.db->sequence_count(), cfg.db->letters() };
	format.print_query_epilog(info);
}

// Called once the last reference to the query state has been released, i.e. all target
// tasks and extensions of the query are done. Completes the output of the query and hands
// it over to the output sink.
void Run::finalize_query(QueryState* state) {
	print_query_epilog(state->query_id, cfg, *output_format, *state->out);
	if (state->hit_num != 0)
		statistics.inc(Statistics::ALIGNED);
	output_sink->push(state->query_id, state->out);
	delete state;

	bool finished;
	{
		std::lock_guard<std::mutex> lock(mtx);
		finished = --queries_pending == 0 && all_queries_submitted;
	}
	if (finished)
		pool.stop();
}

void Run::release(QueryState* state) {
	if (state->refs.fetch_sub(1, std::memory_order_acq_rel) == 1)
		finalize_query(state);
}

void Run::process_extension(Extension& e) {
	extend(e.state, e.target, cfg, output_format.get());
	release(e.state);
}

// Task handling one target of a query, starting at the first hit of the target. Chains the
// task for the next target of the same query before processing its own target. The
// reference to the query state held by this task is handed over to the extension it
// submits, so the task must not touch the state after the submission.
void Run::process_target(QueryState* state, Search::Hit* begin) {
	const SequenceSet& targets = cfg.target->seqs();
	const BlockId target = targets.local_position((uint64_t)begin->subject_).first;
	const uint64_t target_begin = targets.position(target, 0), target_end = targets.position(target + 1, 0);
	Search::Hit* end = begin + 1;
	while (end < state->hits_end && (uint64_t)end->subject_ < target_end)
		++end;

	if (end < state->hits_end) {
		state->refs.fetch_add(1, std::memory_order_relaxed);
		pool.submit([this, state, end] { process_target(state, end); });
	}

	std::vector<::Extension::SeedHit> seed_hits;
	seed_hits.reserve(end - begin);
	for (const Search::Hit* hit = begin; hit < end; ++hit)
		seed_hits.push_back({ (int)hit->seed_offset_, (int)((uint64_t)hit->subject_ - target_begin), hit->score_, hit->frame() });
	Statistics stats;
	::Extension::ungapped_stage(seed_hits.begin(), seed_hits.end(), state->query, target, 0, stats, *cfg.target, cfg.extension_mode,
		*std::pmr::get_default_resource(), cfg);
	statistics += stats;

	pool.submit_extension(Extension(state, target));
}

// Task handling one query, starting at the first hit of the query. Chains the task for the
// next query, then sets up the query state and continues as the task of the first target.
void Run::process_query(Search::Hit* begin) {
	Search::Hit* const hits_end = hits + hit_count;
	const BlockId query = begin->query_;
	Search::Hit* end = begin;
	while (end < hits_end && end->query_ == query)
		++end;

	if (end < hits_end) {
		{
			std::lock_guard<std::mutex> lock(mtx);
			++queries_pending;
		}
		pool.submit([this, end] { process_query(end); });
		for (BlockId query_id = query + 1; query_id < end->query_; ++query_id)
			output_sink->push(query_id, nullptr);
	}
	else {
		std::lock_guard<std::mutex> lock(mtx);
		all_queries_submitted = true;
	}

	std::sort(begin, end, Search::Hit::CmpSubject());
	Statistics stats;
	QueryState* state = new QueryState(this, query, end, stats);
	statistics += stats;
	print_query_intro(query, cfg, *output_format, *state->out);
	process_target(state, begin);
}

void run(Search::Hit* hits, uint64_t hit_count, const Search::Config& cfg) {
	if (hit_count == 0)
		return;

	Run run(hits, hit_count, cfg);

	for (BlockId query_id = static_cast<BlockId>(output_sink->begin()); query_id < hits->query_; ++query_id)
		output_sink->push(query_id, nullptr);

	run.pool.submit([&run] { run.process_query(run.hits); });
}

}

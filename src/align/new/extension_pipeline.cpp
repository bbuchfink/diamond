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

#include <cassert>
#include <cstdint>
#include <limits>
#include <new>
#include <utility>
#include <string>
#include <algorithm>
#include <array>
#include <fstream>
#include <stdexcept>
#include "basic/match.h"
#include "basic/statistics.h"
#include "output/output_format.h"
#include "data/block/block.h"
#include "data/sequence_file.h"
#include "output/output.h"
#include "align/target.h"
#include "align/culling.h"
#include "dp/anchored_swipe/config.h"
#include "thread_pool.h"
#include "extension_pipeline.h"
#include "dp/dp.h"

using std::list;

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

// Conservative upper bound on the band of any extension of a query and a target of the current blocks,
// independent of the anchors. The anchored swipe clips the diagonal range of an extension to the
// sub-matrix it extends into, which is at most as large as the whole matrix of the pair, so the band
// never exceeds the qlen + tlen - 1 diagonals of that matrix.
static Loc block_band_max(const Search::Config& cfg) {
	const SequenceSet& queries = cfg.query->seqs(), & targets = cfg.target->seqs();
	const int64_t query_len = queries.max_len(0, queries.size()), target_len = targets.max_len(0, targets.size());
	return (Loc)std::max(query_len + target_len - 1, (int64_t)0);
}

Run::Run(Search::Hit* hits, uint64_t hit_count, const Search::Config& cfg):
	hits(hits),
	hit_count(hit_count),
	cfg(cfg),
	output_format(cfg.output_format->clone()),
	band_max(block_band_max(cfg)),
	score_tables(::Stats::CBS::matrix_adjust(static_cast<unsigned>(config.comp_based_stats_.get(::Stats::DEFAULT_CBS)))),
	queries_pending(1),
	all_queries_submitted(false),
	output_callback([this](std::vector<Extension>& extensions) { this->output_extensions(extensions); }),
	pool(config.threads_, config.swipe_batch, band_max, [this](std::vector<Extension>& batch, ExtensionQueue& queue) { this->process_extensions(batch, queue); })
{
}

// Prints an alignment of the query of `state` and a target. It carries no traceback, only the
// score, e-value and the query/target ranges.
static void print_alignment(QueryState* state, BlockId target_id, const Sequence& target, const Hsp& hsp,
	const Search::Config& cfg, OutputFormat* format, Statistics& stats)
{
	const BlockId query_id = state->query_id;
	const Sequence& query = state->query.sequence.front();

	const OId target_oid = cfg.target->block_id2oid(target_id);
	const bool all_seqids = flag_any(format->flags, Output::Flags::ALL_SEQIDS);
	const std::string target_title = cfg.target->has_ids()
		? cfg.target->ids()[target_id]
		: (flag_any(format->flags, Output::Flags::SSEQID) ? cfg.db->seqid(target_oid, all_seqids, true) : "");
	const Sequence& output_target = cfg.target->unmasked_seqs().empty() ? target : cfg.target->unmasked_seqs()[target_id];
	const double query_self_score = cfg.query->has_self_aln() ? cfg.query->self_aln_score(query_id) : 0.0;
	const double target_self_score = cfg.target->has_self_aln() ? cfg.target->self_aln_score(target_id) : 0.0;

	std::lock_guard<SpinLock> lock(state->mtx);
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
	stats.inc(Statistics::MATCHES);
	stats.inc(Statistics::PAIRWISE);
}

// Applies the e-value cutoff and the filters of the regular pipeline (coverage, approximate identity, self hits)
// to an alignment of the query of `state` and a target, and prints it if it passes.
static void report(QueryState* state, BlockId target_id, const Sequence& target, Hsp& hsp, const Search::Config& cfg, OutputFormat* format, Statistics& stats) {
	if (hsp.evalue > config.max_evalue)
		return;
	const ::Extension::Query& query = state->query;
	const char* target_title = config.no_self_hits ? cfg.target->ids()[target_id] : nullptr;
	const double target_self_score = cfg.target->has_self_aln() ? cfg.target->self_aln_score(target_id) : 0.0;
	if (::Extension::filter_hsp(hsp, query.source_length, query.title, target.length(), target_title, query.sequence.front(), target,
		query.self_alignment_score, target_self_score, format))
		return;
	print_alignment(state, target_id, target, hsp, cfg, format, stats);
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
	if (state->hit_num != 0) {
		Statistics stats;
		stats.inc(Statistics::ALIGNED);
		statistics += stats;
	}
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

// Widens the diagonal ranges of the anchor for both extension directions by the band on either side,
// so that the anchored swipe computes the extensions in the widened ranges. Unset ranges remain empty.
// The band is the one of the banded extension of the regular pipeline, which depends on the query length.
static Anchor add_band(Anchor anchor, Loc query_len, ::Extension::Mode mode) {
	const Loc band = ::Extension::band(query_len, mode);
	anchor.d_min_left -= band;
	anchor.d_max_left += band;
	anchor.d_min_right -= band;
	anchor.d_max_right += band;
	return anchor;
}

// A direction is only extended if at least MIN_LETTERS letters of both the query and the target
// are left on that side of the anchor.
static bool extend_right(const TargetPair& pair) {
	return std::min(pair.query.length() - pair.anchor.query_end(), pair.target.length() - pair.anchor.subject_end()) >= DpTarget::MIN_LETTERS;
}

static bool extend_left(const TargetPair& pair) {
	return std::min(pair.anchor.query_begin(), pair.anchor.subject_begin()) >= DpTarget::MIN_LETTERS;
}

// Rescores the ungapped anchor under a score table (see TargetPair::score_table) and trims it to its highest
// scoring segment of at least one letter. The chaining found the anchor using the unadjusted matrix, so under the
// adjusted one its ends may score negatively. Trimming does not lose any score, as the extensions of the anchor
// can still include the trimmed letters. The diagonal ranges of the anchor are kept.
static void rescore_anchor(Anchor& anchor, const Sequence& query, const Sequence& target, const int8_t* score_table) {
	Score score = 0, best = std::numeric_limits<Score>::min();
	Loc begin = 0, best_begin = 0, best_end = 0;
	for (Loc k = 0; k < anchor.len; ++k) {
		if (score <= 0) {
			score = 0;
			begin = k;
		}
		score += score_table[((int)letter_mask(target[anchor.j + k]) << 5) + (int)letter_mask(query[anchor.i + k])];
		if (score > best) {
			best = score;
			best_begin = begin;
			best_end = k + 1;
		}
	}
	anchor.i += best_begin;
	anchor.j += best_begin;
	anchor.len = best_end - best_begin;
	anchor.score = best;
}

// Computes a batch of extensions popped from one of the extension queues, possibly belonging to different
// queries, in one call of the anchored swipe, which passes them to the output callback when done. The
// anchored swipe also gets the queue itself, so it keeps popping extensions submitted to it while it is
// running. This is why the band limit does not depend on the extensions of the batch, but is the limit of
// the queue, which covers any extension it may receive.
void Run::process_extensions(std::vector<Extension>& batch, ExtensionQueue& queue) {
	Statistics stats;
	const DP::AnchoredSwipe::Config swipe_cfg{ 0, stats, cfg.extension_mode, queue.band };
	DP::AnchoredSwipe::anchored_swipe(batch, queue, swipe_cfg, output_callback);
	statistics += stats;
}

// Merges the results of the (up to two) extensions of a pair into one alignment, applies the
// filters and prints it. Called once all extensions of the pair have been computed, or directly
// by the target task if the anchor is not extended into either direction. The reference to the
// query state held by the pair is released here.
void Run::finish_pair(TargetPair* pair, Statistics& stats) {
	const Anchor& anchor = pair->anchor;
	// A direction that was not extended contributes zero.
	const Score score = anchor.score + pair->score[RIGHT] + pair->score[LEFT];
	const Interval query_range(anchor.query_begin() - pair->query_len[LEFT], anchor.query_end() + pair->query_len[RIGHT]),
		target_range(anchor.subject_begin() - pair->target_len[LEFT], anchor.subject_end() + pair->target_len[RIGHT]);

	stats.inc(Statistics::EXT16);
	const Loc qlen = pair->query.length(), tlen = pair->target.length();
	Hsp hsp;
	hsp.score = score;
	hsp.evalue = score_matrix.evalue(score, qlen, tlen); // consider adjusted matrix
	hsp.bit_score = score_matrix.bitscore(score);
	hsp.corrected_bit_score = score_matrix.bitscore_corrected(score, qlen, tlen);
	hsp.query_range = hsp.query_source_range = query_range;
	hsp.subject_range = hsp.subject_source_range = target_range;
	hsp.approx_id = hsp.approx_id_percent(pair->query, pair->target);
	report(pair->state, pair->target_id, pair->target, hsp, cfg, output_format.get(), stats);

	QueryState* const state = pair->state;
	PairBlock* const block = pair->block;
	pair->~TargetPair();
	block->release();
	release(state);
}

// Output callback of the anchored swipe. Every extension stores its result in the shared state
// of its pair; the extension of the pair arriving last merges the results and outputs them.
void Run::output_extensions(std::vector<Extension>& extensions) {
	Statistics stats;
	for (const Extension& e : extensions)
		if (e.finish())
			finish_pair(e.pair, stats);
	statistics += stats;
}

// Returns the end of the hits of the target whose first hit is begin.
static Search::Hit* target_end(Search::Hit* begin, Search::Hit* hits_end, const SequenceSet& targets) {
	const BlockId target = targets.local_position((uint64_t)begin->subject_).first;
	const uint64_t end_pos = targets.position(target + 1, 0);
	Search::Hit* end = begin + 1;
	while (end < hits_end && (uint64_t)end->subject_ < end_pos)
		++end;
	return end;
}

// Computes the anchor of one target of a query from its hits [begin, end) and submits the extensions
// of the anchor to the pool. The pair created for the target holds a reference to the query state of its own.
void Run::process_target(QueryState* state, Search::Hit* begin, Search::Hit* end, PairBlock* block, Statistics& stats) {
	const SequenceSet& targets = cfg.target->seqs();
	const BlockId target = targets.local_position((uint64_t)begin->subject_).first;
	const uint64_t target_begin = targets.position(target, 0);

	std::vector<::Extension::SeedHit> seed_hits;
	seed_hits.reserve(end - begin);
	for (const Search::Hit* hit = begin; hit < end; ++hit)
		seed_hits.push_back({ (int)hit->seed_offset_, (int)((uint64_t)hit->subject_ - target_begin), hit->score_, hit->frame() });
	const ::Extension::WorkTarget work_target = ::Extension::ungapped_stage(seed_hits.begin(), seed_hits.end(), state->query, target, 0, stats, *cfg.target,
		cfg.extension_mode, *std::pmr::get_default_resource(), cfg);
	const list<ApproxHsp>& hsps = work_target.hsp[0];
	if (work_target.done) {
		// The ungapped stage already found an alignment passing the filters (--hamming-ext), which is final
		// like in the regular pipeline. Its anchor carries no diagonal ranges to extend.
		Hsp hsp(hsps.front(), state->query.sequence.front().length(), targets[target].length());
		report(state, target, targets[target], hsp, cfg, output_format.get(), stats);
		return;
	}
	const Sequence& query = state->query.sequence.front();
	Anchor anchor = hsps.empty() ? Anchor() : add_band(hsps.front().max_diag, query.length(), cfg.extension_mode);
	if (anchor.score <= 0)
		return;

	// The ungapped stage has computed the compositionally adjusted matrix of the pair, if the composition based
	// statistics call for one. The extensions are computed with it, and the anchor is rescored with it, so that the
	// score of the alignment is the one under the adjusted matrix.
	const int8_t* score_table = nullptr;
	if (work_target.matrix) {
		int8_t* table = block->next_score_table();
		DP::AnchoredSwipe::make_score_table(work_target.matrix->scores.data(), AMINO_ACID_COUNT, table);
		rescore_anchor(anchor, query, targets[target], table);
		score_table = table;
	}

	// The two directions of the anchor are extended independently of each other, so they are
	// submitted as two extensions sharing the state of the pair.
	state->refs.fetch_add(1, std::memory_order_relaxed);
	TargetPair* pair = block->new_pair(state, target, query, targets[target], anchor, score_table);
	const bool right = extend_right(*pair), left = extend_left(*pair);
	// Set before the submission: an extension may already be computed and finish the pair
	// while the other one is still being submitted. The pair must not be accessed after the last submission.
	pair->pending.store(int(right) + int(left), std::memory_order_relaxed);
	if (right)
		pool.submit_extension(Extension(pair, RIGHT, anchor.query_end(), anchor.subject_end(), anchor.d_min_right, anchor.d_max_right));
	if (left)
		pool.submit_extension(Extension(pair, LEFT, anchor.query_begin() - 1, anchor.subject_begin() - 1, anchor.d_min_left, anchor.d_max_left));
	if (!right && !left)
		// Nothing to extend, the anchor itself is the alignment.
		finish_pair(pair, stats);
}

// Task handling a chunk of up to TARGETS_PER_TASK targets of a query, starting at the first hit of a target.
// Chains the task for the next chunk of the same query before processing its own targets, so that the
// targets of one query are processed in parallel. Every extension is submitted to the pool as soon as it has
// been created. The task holds a reference to the query state, which it releases when done.
void Run::process_targets(QueryState* state, Search::Hit* begin) {
	const SequenceSet& targets = cfg.target->seqs();
	std::array<Search::Hit*, TARGETS_PER_TASK + 1> bounds;
	bounds[0] = begin;
	int n = 0;
	while (n < TARGETS_PER_TASK && bounds[n] < state->hits_end) {
		bounds[n + 1] = target_end(bounds[n], state->hits_end, targets);
		++n;
	}

	if (bounds[n] < state->hits_end) {
		state->refs.fetch_add(1, std::memory_order_relaxed);
		Search::Hit* const next = bounds[n];
		pool.submit([this, state, next] { process_targets(state, next); });
	}

	Statistics stats;
	PairBlock* const block = PairBlock::create(n, score_tables);
	for (int k = 0; k < n; ++k)
		process_target(state, bounds[k], bounds[k + 1], block, stats);
	statistics += stats;
	block->release();
	release(state);
}

// Task handling one query, starting at the first hit of the query. Chains the task for the
// next query, then sets up the query state and continues as the task of the first chunk of targets.
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
	process_targets(state, begin);
}

static std::mutex band_histogram_mtx;
static std::vector<int64_t> band_histogram;

void add_band_histogram(const std::vector<int64_t>& counts) {
	std::lock_guard<std::mutex> lock(band_histogram_mtx);
	if (band_histogram.size() < counts.size())
		band_histogram.resize(counts.size(), 0);
	for (size_t b = 0; b < counts.size(); ++b)
		band_histogram[b] += counts[b];
}

void write_band_histogram(const std::string& file) {
	std::ofstream out(file);
	if (!out)
		throw std::runtime_error("Error opening file " + file);
	std::lock_guard<std::mutex> lock(band_histogram_mtx);
	for (size_t b = 0; b < band_histogram.size(); ++b)
		if (band_histogram[b] != 0)
			out << b << '\t' << band_histogram[b] << '\n';
	if (!out)
		throw std::runtime_error("Error writing file " + file);
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

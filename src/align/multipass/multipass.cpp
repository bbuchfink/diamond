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

#include <algorithm>
#include <array>
#include <cassert>
#include <climits>
#include <limits>
#include <list>
#include <vector>
#include <string>
#include "basic/config.h"
#include "basic/value.h"
#include "basic/match.h"
#include "data/block/block.h"
#include "data/sequence_set.h"
#include "data/sequence_file.h"
#include "output/output.h"
#include "output/output_format.h"
#include "stats/score_matrix.h"
#include "align/culling.h"
#include "util/memory/memory_resource.h"
#include "align/target.h"
#include "align/def.h"
#include "dp/dp.h"
#include "dp/anchored_swipe/config.h"
#include "anchored_swipe.h"
#include "multipass.h"

namespace MultipassPipeline {

// The hits carry the id of the query context (frame), all contexts of a query belong to the same query task.
static BlockId query_id(const Search::Hit& hit) {
	return hit.query_ / align_mode.query_contexts;
}

QueryState::QueryState(BlockId query_id, Search::Hit* hits_end, Statistics& stats, const Search::Config& cfg) :
	query_id(query_id),
	hits_end(hits_end),
	query(query_id, stats, cfg, *std::pmr::get_default_resource()),
	refs(1)
{
}

Run::Run(Search::Hit* hits, uint64_t hit_count, const Search::Config& cfg, File* output_file) :
	hits(hits),
	hit_count(hit_count),
	cfg(cfg),
	output_format(cfg.output_format->clone()),
	output_file(output_file),
	query_aligned(new std::atomic<bool>[cfg.query->seqs().size() / align_mode.query_contexts]()),
	queries_pending(1),
	all_queries_submitted(false),
	extension_buffers(std::max(config.threads_, 1)),
	extension_count(0),
	pool(config.threads_)
{
	// Stats::TargetMatrix takes the scores of the letters that it does not adjust from the unadjusted matrix.
	std::array<int8_t, 32 * AMINO_ACID_COUNT> matrix;
	matrix.fill(0);
	for (size_t i = 0; i < AMINO_ACID_COUNT; ++i)
		for (size_t j = 0; j < AMINO_ACID_COUNT; ++j)
			matrix[i * 32 + j] = (int8_t)std::clamp(score_matrix(i, j) * config.cbs_matrix_scale, SCHAR_MIN, SCHAR_MAX);
	DP::AnchoredSwipe::make_score_table(matrix.data(), AMINO_ACID_COUNT, common_scores);
}

// Called once the last reference to the query state has been released, i.e. all target
// tasks of the query are done. Stops the pool if this was the last query.
void Run::finalize_query(QueryState* state) {
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

// Returns the end of the hits of the target whose first hit is begin.
static Search::Hit* target_end(Search::Hit* begin, Search::Hit* hits_end, const SequenceSet& targets) {
	const BlockId target = (BlockId)targets.local_position((uint64_t)begin->subject_).first;
	const uint64_t end_pos = targets.position(target + 1, 0);
	Search::Hit* end = begin + 1;
	while (end < hits_end && (uint64_t)end->subject_ < end_pos)
		++end;
	return end;
}

// Widens the diagonal ranges of the anchor for both extension directions by the band on either side,
// so that the extensions are computed in the widened ranges. Unset ranges remain empty.
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

// Appends the formatted output of a task to the output file and clears it.
void Run::write(TextBuffer& out) {
	if (out.size() == 0)
		return;
	{
		std::lock_guard<std::mutex> lock(output_mtx);
		output_file->write(out.data(), out.size());
	}
	out.clear();
}

// Applies the e-value cutoff and the filters of the regular pipeline (coverage, approximate identity, self hits)
// to an alignment of a query and a target, and prints it to out if it passes. The alignment carries no traceback,
// only the score, e-value and the query/target ranges.
void Run::report(BlockId query_id, const Sequence& query, BlockId target_id, const Sequence& target, Hsp& hsp, TextBuffer& out, Statistics& stats) {
	if (hsp.evalue > config.max_evalue)
		return;
	const char* const query_title = cfg.query->ids()[query_id];
	const int query_source_len = align_mode.query_translated ? (int)cfg.query->source_seqs()[query_id].length() : (int)query.length();
	const double query_self_score = cfg.query->has_self_aln() ? cfg.query->self_aln_score(query_id) : 0.0;
	const double target_self_score = cfg.target->has_self_aln() ? cfg.target->self_aln_score(target_id) : 0.0;
	if (::Extension::filter_hsp(hsp, query_source_len, query_title, target.length(), config.no_self_hits ? cfg.target->ids()[target_id] : nullptr,
		query, target, query_self_score, target_self_score, output_format.get()))
		return;

	const OId target_oid = cfg.target->block_id2oid(target_id);
	const bool all_seqids = flag_any(output_format->flags, Output::Flags::ALL_SEQIDS);
	const std::string target_title = cfg.target->has_ids()
		? cfg.target->ids()[target_id]
		: (flag_any(output_format->flags, Output::Flags::SSEQID) ? cfg.db->seqid(target_oid, all_seqids, true) : "");
	const Sequence& output_target = cfg.target->unmasked_seqs().empty() ? target : cfg.target->unmasked_seqs()[target_id];
	Output::Info info{ cfg.query->seq_info(query_id), false, cfg.db.get(), out, Util::Seq::AccessionParsing(), cfg.db->sequence_count(), cfg.db->letters() };
	// The alignments of a query are not numbered, as they are not collected per query.
	output_format->print_match(HspContext(hsp,
		query_id,
		cfg.query->block_id2oid(query_id),
		TranslatedSequence(query),
		query_title,
		target_oid,
		static_cast<unsigned>(target.length()),
		target_title.c_str(),
		0,
		0,
		output_target,
		0,
		query_self_score,
		target_self_score), info);
	stats.inc(Statistics::MATCHES);
	stats.inc(Statistics::PAIRWISE);
	// Read first, as the flag of an aligned query is not written again, which keeps its cache line shared.
	std::atomic<bool>& aligned = query_aligned[query_id];
	if (!aligned.load(std::memory_order_relaxed) && !aligned.exchange(true, std::memory_order_relaxed))
		stats.inc(Statistics::ALIGNED);
}

// Merges the results of the (up to two) extensions of a pair into one alignment, applies the filters and prints
// it to out. Called once all extensions of the pair have been computed, or directly by the target task if the
// anchor is not extended into either direction.
void Run::finish_pair(const TargetPair& pair, TextBuffer& out, Statistics& stats) {
	const Anchor& anchor = pair.anchor;
	// A direction that was not extended contributes zero.
	const Score score = anchor.score + pair.score[RIGHT] + pair.score[LEFT];
	const Interval query_range(anchor.query_begin() - pair.query_len[LEFT], anchor.query_end() + pair.query_len[RIGHT]),
		target_range(anchor.subject_begin() - pair.target_len[LEFT], anchor.subject_end() + pair.target_len[RIGHT]);

	stats.inc(Statistics::EXT16);
	const Loc qlen = pair.query.length(), tlen = pair.target.length();
	Hsp hsp;
	hsp.score = score;
	hsp.evalue = score_matrix.evalue(score, qlen, tlen);
	hsp.bit_score = score_matrix.bitscore(score);
	hsp.corrected_bit_score = score_matrix.bitscore_corrected(score, qlen, tlen);
	hsp.query_range = hsp.query_source_range = query_range;
	hsp.subject_range = hsp.subject_source_range = target_range;
	hsp.approx_id = hsp.approx_id_percent(pair.query, pair.target);
	report(pair.query_id, pair.query, pair.target_id, pair.target, hsp, out, stats);
}

// Computes the anchor of one target of a query from its hits [begin, end) and appends the pair and the
// extensions of the anchor to the buffer of the calling worker, where they stay until the end of the run.
// A pair that has nothing to extend is printed to out right away.
void Run::process_target(QueryState* state, Search::Hit* begin, Search::Hit* end, ExtensionBuffer& buffer, TextBuffer& out, Statistics& stats) {
	const SequenceSet& targets = cfg.target->seqs();
	const BlockId target = (BlockId)targets.local_position((uint64_t)begin->subject_).first;
	const uint64_t target_begin = targets.position(target, 0);

	std::vector<::Extension::SeedHit> seed_hits;
	seed_hits.reserve(end - begin);
	for (const Search::Hit* hit = begin; hit < end; ++hit)
		seed_hits.push_back({ (int)hit->seed_offset_, (int)((uint64_t)hit->subject_ - target_begin), hit->score_, hit->frame() });
	const ::Extension::WorkTarget work_target = ::Extension::ungapped_stage(seed_hits.begin(), seed_hits.end(), state->query, target, 0, stats, *cfg.target,
		cfg.extension_mode, *std::pmr::get_default_resource(), cfg);
	const std::list<ApproxHsp>& hsps = work_target.hsp[0];
	const Sequence& query = state->query.sequence.front();
	if (work_target.done) {
		// The ungapped stage already found an alignment passing the filters (--hamming-ext), which is final
		// like in the regular pipeline. There is nothing to extend.
		Hsp hsp(hsps.front(), query.length(), targets[target].length());
		report(state->query_id, query, target, targets[target], hsp, out, stats);
		return;
	}
	Anchor anchor = hsps.empty() ? Anchor() : add_band(hsps.front().max_diag, query.length(), cfg.extension_mode);
	if (anchor.score <= 0)
		return;

	// The ungapped stage has computed the compositionally adjusted matrix of the pair, if the composition based
	// statistics call for one. The anchor is rescored with it, so that its score is the one under the adjusted matrix.
	std::pair<const int8_t*, bool> score_table(nullptr, false);
	if (work_target.matrix) {
		int8_t table[32 * 32];
		DP::AnchoredSwipe::make_score_table(work_target.matrix->scores.data(), AMINO_ACID_COUNT, table);
		rescore_anchor(anchor, query, targets[target], table);
		score_table = buffer.add_score_table(table, common_scores);
	}

	// The two directions of the anchor are extended independently of each other, so they are buffered as two
	// extensions referring to the pair, which is printed once both have been computed.
	TargetPair& pair = buffer.pairs.emplace_back(state->query_id, target, query, targets[target], anchor, score_table.first, score_table.second);
	const bool right = extend_right(pair), left = extend_left(pair);
	pair.pending.store(int(right) + int(left), std::memory_order_relaxed);
	if (right)
		buffer.add_extension(&pair, RIGHT, anchor.query_end(), anchor.subject_end(), anchor.d_min_right, anchor.d_max_right);
	if (left)
		buffer.add_extension(&pair, LEFT, anchor.query_begin() - 1, anchor.subject_begin() - 1, anchor.d_min_left, anchor.d_max_left);
	if (!right && !left)
		// Nothing to extend, the anchor itself is the alignment.
		finish_pair(pair, out, stats);
}

// Task handling a batch of up to TARGETS_PER_TASK targets of a query, starting at the first hit of a target.
// Chains the task for the next batch of the same query before processing its own targets, so that the
// targets of one query are processed in parallel. The task holds a reference to the query state, which it
// releases when done.
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
	assert(ThreadPool::worker_index() >= 0 && (size_t)ThreadPool::worker_index() < extension_buffers.size());
	ExtensionBuffer& buffer = extension_buffers[ThreadPool::worker_index()];
	TextBuffer out;
	for (int k = 0; k < n; ++k)
		process_target(state, bounds[k], bounds[k + 1], buffer, out, stats);
	write(out);
	statistics += stats;
	release(state);
}

// Task handling one query, starting at the first hit of the query. Chains the task for the
// next query, then sets up the query state and continues as the task of the first batch of targets.
void Run::process_query(Search::Hit* begin) {
	Search::Hit* const hits_end = hits + hit_count;
	const BlockId query = query_id(*begin);
	Search::Hit* end = begin;
	while (end < hits_end && query_id(*end) == query)
		++end;

	if (end < hits_end) {
		{
			std::lock_guard<std::mutex> lock(mtx);
			++queries_pending;
		}
		pool.submit([this, end] { process_query(end); });
	}
	else {
		std::lock_guard<std::mutex> lock(mtx);
		all_queries_submitted = true;
	}

	std::sort(begin, end, Search::Hit::CmpSubject());
	Statistics stats;
	QueryState* state = new QueryState(query, end, stats, cfg);
	statistics += stats;
	process_targets(state, begin);
}

// Orders the buffered extensions by band. Requires the workers of the pool to have joined.
// This is a counting sort on the bands, which the buffers have counted the extensions of. As the position of
// every extension in the order follows from the counts, each buffer is distributed by a task of its own. Only
// pointers to the extensions are moved.
void Run::sort_extensions() {
	size_t bands = 0;
	for (const ExtensionBuffer& buffer : extension_buffers)
		bands = std::max(bands, buffer.band_counts.size());

	// Element [w][b] is the position of the next extension of band b of buffer w.
	std::vector<std::vector<size_t>> positions(extension_buffers.size(), std::vector<size_t>(bands));
	size_t n = 0;
	for (size_t b = 0; b < bands; ++b)
		for (size_t w = 0; w < extension_buffers.size(); ++w) {
			positions[w][b] = n;
			if (b < extension_buffers[w].band_counts.size())
				n += extension_buffers[w].band_counts[b];
		}
	extension_count = n;
	sorted_extensions.reset(new const Extension*[n]);

	ThreadPool sort_pool(config.threads_);
	for (size_t w = 0; w < extension_buffers.size(); ++w)
		sort_pool.submit([this, w, &positions] {
			size_t* const pos = positions[w].data();
			const Extension** const sorted = sorted_extensions.get();
			extension_buffers[w].extensions.for_each_chunk([pos, sorted](const Extension* begin, const Extension* end) {
				for (const Extension* e = begin; e < end; ++e)
					sorted[pos[e->band]++] = e;
			});
		});
	sort_pool.stop();
	sort_pool.join();
}

// Computes the buffered extensions with the anchored swipe. Requires the extensions to be sorted by band.
// They are split into packages of PACKAGE_SIZE extensions that are consecutive in that order, i.e. have similar
// bands. Each package is computed by one task in one call of the kernel, with the widest band of the package as
// the band limit, so that the kernel neither computes the extensions over a band much wider than their own nor
// sizes its matrix and the query paddings for the widest band of the run. Every extension stores its result in its pair. The extension of a pair computed last merges the
// results and prints the alignment, which may happen in a different task than for the other extension of the pair.
void Run::compute_extensions() {
	ThreadPool swipe_pool(config.threads_);
	const Extension* const* const sorted = sorted_extensions.get();
	// The packages are queued in the order of descending band. The cost of a package grows with its band, and
	// starting with the expensive ones lets the cheap ones even out the load of the workers at the end.
	for (size_t end = extension_count; end > 0; end -= std::min(PACKAGE_SIZE, end)) {
		const size_t begin = end - std::min(PACKAGE_SIZE, end);
		swipe_pool.submit([this, sorted, begin, end] {
			Statistics stats;
			TextBuffer out;
			// The last extension of the package has the widest band.
			const Loc band_max = sorted[end - 1]->band;
			anchored_swipe(sorted + begin, sorted + end, band_max, common_scores, stats, [this, &out, &stats](const Extension& e, Score score, Loc query_len, Loc target_len) {
				if (!e.pair->finish(e.direction, score, query_len, target_len))
					return;
				finish_pair(*e.pair, out, stats);
				if (out.size() >= OUTPUT_FLUSH_SIZE)
					write(out);
			});
			write(out);
			statistics += stats;
		});
	}
	// The workers complete the tasks queued before they exit, and the pool waits for them in its destructor.
	swipe_pool.stop();
}

void run(Search::Hit* hits, uint64_t hit_count, const Search::Config& cfg, File* output_file) {
	if (hit_count == 0)
		return;

	Run run(hits, hit_count, cfg, output_file);
	{
		TaskTimer timer("Computing chaining");
		run.pool.submit([&run] { run.process_query(run.hits); });
		run.pool.join();
	}
	{
		TaskTimer timer("Sorting extensions");
		run.sort_extensions();
	}
	{
		TaskTimer timer("Computing extensions");
		run.compute_extensions();
	}
}

}

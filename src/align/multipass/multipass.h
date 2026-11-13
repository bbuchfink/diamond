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

#pragma once
#include <atomic>
#include <cstdint>
#include <memory>
#include <mutex>
#include <vector>
#include "search/hit.h"
#include "run/config.h"
#include "basic/statistics.h"
#include "basic/match.h"
#include "util/io/file.h"
#include "util/text_buffer.h"
#include "output/output_format.h"
#include "align/extend.h"
#include "extension.h"
#include "thread_pool.h"

namespace MultipassPipeline {

// Heap allocated state shared by the target tasks of one query, each of which holds a reference
// to it. It is finalized and freed once the last reference has been released.
struct QueryState {

	QueryState(BlockId query_id, Search::Hit* hits_end, Statistics& stats, const Search::Config& cfg);

	const BlockId query_id;
	Search::Hit* const hits_end;  // end of the hits of this query (sorted by subject)
	const ::Extension::Query query;
	std::atomic<int64_t> refs;    // target tasks in flight
};

struct Run {

	// Maximum number of targets of a query handled by one target task.
	static constexpr int TARGETS_PER_TASK = 64;

	// Number of extensions computed in one call of the anchored swipe kernel (see compute_extensions). Smaller
	// packages even out the load of the workers better, as the cost of a package varies a lot with its band and
	// the lengths of its sequences. Larger ones leave less of the kernel idle while its last extensions complete.
	static constexpr size_t PACKAGE_SIZE = 1024;

	// Size of the formatted output of a task at which it is written to the output file.
	static constexpr size_t OUTPUT_FLUSH_SIZE = (size_t)1 << 20;

	Run(Search::Hit* hits, uint64_t hit_count, const Search::Config& cfg, File* output_file);
	void process_query(Search::Hit* begin);
	void process_targets(QueryState* state, Search::Hit* begin);
	void process_target(QueryState* state, Search::Hit* begin, Search::Hit* end, ExtensionBuffer& buffer, TextBuffer& out, Statistics& stats);
	void release(QueryState* state);
	void finalize_query(QueryState* state);
	void sort_extensions();
	void compute_extensions();
	void finish_pair(const TargetPair& pair, TextBuffer& out, Statistics& stats);
	void report(BlockId query_id, const Sequence& query, BlockId target_id, const Sequence& target, Hsp& hsp, TextBuffer& out, Statistics& stats);
	void write(TextBuffer& out);

	Search::Hit* hits;
	uint64_t hit_count;
	const Search::Config& cfg;
	std::unique_ptr<OutputFormat> output_format;

	// The tasks format their alignments into buffers of their own and append them to the output file, in no
	// particular order and with the alignments of different queries mixed.
	File* const output_file;
	std::mutex output_mtx;        // guards output_file
	// Set for the queries of the block that an alignment has been printed for, to count the aligned queries.
	const std::unique_ptr<std::atomic<bool>[]> query_aligned;

	std::mutex mtx;               // guards the two members below
	int64_t queries_pending;      // query tasks submitted but whose query is not yet finalized
	bool all_queries_submitted;

	// The scores of the score table of a pair (32x32, see TargetPair::score_table) that are the same for all pairs
	// with a compositionally adjusted matrix, i.e. the ones that a CompactScoreTable does not hold.
	alignas(32) int8_t common_scores[32 * 32];

	// The extensions created over the whole run, indexed by the worker that created them. Computed once the
	// workers have joined (see compute_extensions), and freed after that, as the pool is destroyed first.
	std::vector<ExtensionBuffer> extension_buffers;
	// Pointers to all extensions of the buffers, ordered by band (see sort_extensions).
	std::unique_ptr<const Extension*[]> sorted_extensions;
	size_t extension_count;

	// Runs the query and target tasks. Has to be the last member (joins the workers in its destructor).
	ThreadPool pool;

};

// Computes the alignments for the hits and writes them to the output file. The alignments are written as they are
// completed, i.e. in no particular order, also not grouped by query.
void run(Search::Hit* hits, uint64_t hit_count, const Search::Config& cfg, File* output_file);

}

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
#include <cassert>
#include <new>
#include <mutex>
#include <memory>
#include <vector>
#include "search/hit.h"
#include "run/config.h"
#include "util/text_buffer.h"
#include "util/parallel/mutex.h"
#include "align/extend.h"
#include "extension.h"
#include "thread_pool.h"

namespace ExtensionPipeline {

struct Run;

// Storage of the pairs created by one target task, allocated in one piece with room for the maximum number of
// pairs of the task. Every pair holds a reference to it, as does the task while it creates the pairs. It is freed
// once the last reference has been released, so the memory is recycled by the allocator instead of every pair of
// the run touching fresh pages. If created with score tables, it also holds a score table (see TargetPair::score_table)
// for every pair, placed behind the pairs.
struct alignas(TargetPair) PairBlock {

	static constexpr size_t SCORE_TABLE_SIZE = 32 * 32;

	static PairBlock* create(int capacity, bool score_tables) {
		const size_t size = sizeof(PairBlock) + capacity * (sizeof(TargetPair) + (score_tables ? SCORE_TABLE_SIZE : 0));
		return new (::operator new(size)) PairBlock(capacity);
	}

	// Requires the block to be created with score tables. Returns the score table of the pair created next.
	int8_t* next_score_table() {
		assert(used < capacity);
		return reinterpret_cast<int8_t*>(reinterpret_cast<TargetPair*>(this + 1) + capacity) + used * SCORE_TABLE_SIZE;
	}

	template<typename... Args>
	TargetPair* new_pair(Args&&... args) {
		assert(used < capacity);
		refs.fetch_add(1, std::memory_order_relaxed);
		return new (reinterpret_cast<TargetPair*>(this + 1) + used++) TargetPair(this, std::forward<Args>(args)...);
	}

	void release() {
		if (refs.fetch_sub(1, std::memory_order_acq_rel) == 1) {
			this->~PairBlock();
			::operator delete(this);
		}
	}

private:

	PairBlock(int capacity) :
		refs(1),
		capacity(capacity),
		used(0)
	{}

	std::atomic<int> refs;        // the task creating the pairs and the pairs not yet finished
	const int capacity;
	int used;                     // slots handed out, only accessed by the creating task

};

// Heap allocated state shared by the target tasks and extensions of one query. Every
// target task and every pair whose extensions are in flight holds a reference to it. It is finalized and
// freed once the last reference has been released.
struct QueryState {

	QueryState(Run* run, BlockId query_id, Search::Hit* hits_end, Statistics& stats);

	Run* const run;
	const BlockId query_id;
	Search::Hit* const hits_end;  // end of the hits of this query (sorted by subject)
	const ::Extension::Query query;
	std::atomic<int64_t> refs;    // target tasks and pairs in flight
	unsigned hit_num;             // number of alignments written so far
	// Guards hit_num and out. Only held to print one alignment, but by all the anchored swipes computing extensions
	// of the query, which would mostly be suspended waiting for a mutex.
	SpinLock mtx;
	TextBuffer* out;              // receives the output of the computed extensions
};

struct Run {

	// Maximum number of targets of a query handled by one target task.
	static constexpr int TARGETS_PER_TASK = 64;

	Run(Search::Hit* hits, uint64_t hit_count, const Search::Config& cfg);
	void process_query(Search::Hit* begin);
	void process_targets(QueryState* state, Search::Hit* begin);
	void process_target(QueryState* state, Search::Hit* begin, Search::Hit* end, PairBlock* block, Statistics& stats);
	void process_extensions(std::vector<Extension>& batch, ExtensionQueue& queue);
	void output_extensions(std::vector<Extension>& extensions);
	void finish_pair(TargetPair* pair, Statistics& stats);
	void release(QueryState* state);
	void finalize_query(QueryState* state);

	Search::Hit* hits;
	uint64_t hit_count;
	const Search::Config& cfg;
	std::unique_ptr<OutputFormat> output_format;
	// Bounds the band of any extension of the current blocks. Band limit of the last extension queue, caps the others.
	const Loc band_max;
	// The composition based statistics adjust the matrix per pair, so the pair blocks hold score tables.
	const bool score_tables;

	std::mutex mtx;               // guards the two members below
	int64_t queries_pending;      // query tasks submitted but whose query is not yet finalized
	bool all_queries_submitted;

	const ExtensionCallback output_callback;
	ThreadPool pool;              // has to be the last member (joins the workers in its destructor)

};

void run(Search::Hit* hits, uint64_t hit_count, const Search::Config& cfg);

}

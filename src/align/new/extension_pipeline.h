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
#include <mutex>
#include <memory>
#include "search/hit.h"
#include "run/config.h"
#include "util/text_buffer.h"
#include "align/extend.h"
#include "thread_pool.h"

namespace ExtensionPipeline {

struct Run;

// Heap allocated state shared by the target tasks and extensions of one query. Every
// target task and extension in flight holds a reference to it. It is finalized and
// freed once the last reference has been released.
struct QueryState {

	QueryState(Run* run, BlockId query_id, Search::Hit* hits_end, Statistics& stats);

	Run* const run;
	const BlockId query_id;
	Search::Hit* const hits_end;  // end of the hits of this query (sorted by subject)
	const ::Extension::Query query;
	std::atomic<int64_t> refs;    // target tasks and extensions in flight
	unsigned hit_num;             // number of alignments written so far
	std::mutex mtx;               // guards hit_num and out
	TextBuffer* out;              // receives the output of the computed extensions
};

struct Extension {
	Extension() :
		state(nullptr),
		target(0)
	{}
	Extension(QueryState* state, BlockId target) :
		state(state),
		target(target)
	{}
	QueryState* state;
	BlockId target;
};

struct Run {

	Run(Search::Hit* hits, uint64_t hit_count, const Search::Config& cfg);
	void process_query(Search::Hit* begin);
	void process_target(QueryState* state, Search::Hit* begin);
	void process_extension(Extension& e);
	void release(QueryState* state);
	void finalize_query(QueryState* state);

	Search::Hit* hits;
	uint64_t hit_count;
	const Search::Config& cfg;
	std::unique_ptr<OutputFormat> output_format;

	std::mutex mtx;               // guards the two members below
	int64_t queries_pending;      // query tasks submitted but whose query is not yet finalized
	bool all_queries_submitted;

	ThreadPool<Extension> pool;   // has to be the last member (joins the workers in its destructor)

};

void run(Search::Hit* hits, uint64_t hit_count, const Search::Config& cfg);

}

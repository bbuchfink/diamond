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
#include <algorithm>
#include <atomic>
#include <cfloat>
#include <functional>
#include <iterator>
#include <memory>
#include <mutex>
#include <string>
#include <vector>
#include "basic/sequence.h"
#include "util/geo/interval.h"
#include "util/hsp/approx_hsp.h"
#include "util/parallel/mutex.h"

namespace ExtensionPipeline {

struct QueryState;
struct PairBlock;

// The two directions an anchor is extended into. Used to index the result slots of TargetPair.
enum Direction { RIGHT = 0, LEFT = 1, DIRECTIONS = 2 };

// Shared state of one query/target pair whose anchor is extended. The two directions are
// submitted to the extension queue as independent extensions, each of which holds a reference
// to the pair. The extension arriving last merges the two results into one alignment.
struct TargetPair {
	TargetPair(PairBlock* block, QueryState* state, BlockId target_id, Sequence query, Sequence target, const Anchor& anchor, const int8_t* score_table) :
		block(block),
		state(state),
		target_id(target_id),
		query(query),
		target(target),
		anchor(anchor),
		score_table(score_table),
		pending(0),
		score{ 0, 0 },
		query_len{ 0, 0 },
		target_len{ 0, 0 }
	{}
	PairBlock* const block;         // storage holding the pair
	QueryState* const state;
	const BlockId target_id;
	const Sequence query, target;
	const Anchor anchor;            // its score is the one under the score table of the pair
	// Compositionally adjusted matrix of the pair in the layout of DP::AnchoredSwipe::Options::score_table, which the
	// extensions are computed with. Held by the block of the pair. nullptr if the matrix is not adjusted.
	const int8_t* const score_table;
	std::atomic<int> pending;       // extensions submitted and not yet computed
	// Results of the extensions, indexed by direction. A direction that is not extended stays zero.
	Score score[DIRECTIONS];
	Loc query_len[DIRECTIONS], target_len[DIRECTIONS];
};

// Extension of an anchor into one direction of a query/target pair, computed by the anchored
// swipe as an independent local alignment. The target task submitting it sets the input
// members, the anchored swipe sets the result members.
struct Extension {
	Extension() :
		pair(nullptr),
		i(0),
		j(0),
		d_min(0),
		d_max(0),
		direction(RIGHT),
		score(0),
		query_len(0),
		target_len(0)
	{}
	Extension(TargetPair* pair, Direction direction, Loc i, Loc j, Loc d_min, Loc d_max) :
		pair(pair),
		query(pair->query),
		target(pair->target),
		i(i),
		j(j),
		d_min(d_min),
		d_max(d_max),
		direction(direction),
		score(0),
		query_len(0),
		target_len(0)
	{}
	// Diagonal of the anchor, which the extension starts out on.
	Loc anchor_diag() const {
		return i - j;
	}
	// Width of the diagonal range of the extension, including the anchor diagonal itself.
	// An unset range counts as empty. Bounds the band the anchored swipe computes the extension in, which
	// is the range clipped to the sub-matrix, and selects the extension queue (see queue_index).
	Loc diag_spread() const {
		const Loc d = anchor_diag();
		return std::max(d_max, d) - std::min(d_min, d);
	}
	// Stores the result in the shared state of the pair and returns true if this was the last
	// extension of the pair, i.e. the caller has to merge and output it.
	bool finish() const {
		pair->score[direction] = score;
		pair->query_len[direction] = query_len;
		pair->target_len[direction] = target_len;
		return pair->pending.fetch_sub(1, std::memory_order_acq_rel) == 1;
	}
	TargetPair* pair;
	Sequence query, target;         // full sequences of the pair
	Loc i, j;                       // first cell: end of the anchor (RIGHT), last letter before it (LEFT)
	Loc d_min, d_max;               // diagonal range to compute
	Direction direction;
	Score score;                    // result: score of the extension, 0 if none was found
	Loc query_len, target_len;      // result: extent of the extension in the query/target
};

// The extensions are partitioned into QUEUE_COUNT queues by their band (Extension::diag_spread()),
// so that each call of the anchored swipe runs with a kernel band fitting the extensions it computes.
// QUEUE_BAND_LIMITS holds the upper limits of the band of the extensions of all queues but the last,
// in ascending order. The last queue takes the extensions exceeding all of them, i.e. any band.
static constexpr Loc QUEUE_BAND_LIMITS[] = { 64, 128, 256, 512 };
static constexpr int QUEUE_COUNT = int(std::size(QUEUE_BAND_LIMITS)) + 1;

// Returns the index of the queue for extensions of the given band: the first queue whose limit covers it.
inline int queue_index(Loc band) {
	int k = 0;
	while (k < QUEUE_COUNT - 1 && band > QUEUE_BAND_LIMITS[k])
		++k;
	return k;
}

// Returns the band limit passed to the anchored swipe for the extensions of queue k. band_max bounds the
// band of any extension of the current blocks, which makes it the limit of the last queue and caps the others.
inline Loc queue_band(int k, Loc band_max) {
	return k < QUEUE_COUNT - 1 ? std::min(QUEUE_BAND_LIMITS[k], band_max) : band_max;
}

// Queue of the extensions of one band range waiting to be computed. The target tasks push every extension
// on its own as soon as they have created it, so it can be popped right away. The thread pool pops a batch of extensions
// to start the anchored swipe with (see ThreadPool). Whenever a lane of its kernel completes after the batch is used up,
// the anchored swipe pops the next extension from the same queue to load the lane with.
// As every extension is pushed and popped on its own, by all threads of the pool, the queue is partitioned into
// shards, each of which is a FIFO queue with a lock of its own. A thread pushes to the shard assigned to it (see
// thread_shard) and pops from the shard it popped from last, moving on to the next non-empty one once that is empty.
// This spreads the threads over the shards, so that they rarely wait for each other and, more importantly, do not
// all write to the same cache lines, which is what limits the throughput of a queue shared by many threads no matter
// how it is synchronized. The price is that the extensions are only popped in the order they were pushed in per shard.
// The extensions of a shard are held in a ring buffer, so that pushing and popping neither allocates nor frees memory
// under the lock, except for the rare growth of the buffer.
struct ExtensionQueue {
	ExtensionQueue(Loc band, int shards) :
		band(band),
		shard_count_(std::max(shards, 1)),
		shards_(new Shard[shard_count_])
	{}
	// The shard assigned to the calling thread, which it pushes to and starts popping from. The value is taken modulo
	// the number of shards of the queue. Threads sharing a shard are safe, but contend for it.
	static int& thread_shard() {
		static thread_local int shard = 0;
		return shard;
	}
	// Appends the extension to the shard of the calling thread.
	void push(Extension&& extension) {
		Shard& shard = shards_[thread_shard() % shard_count_];
		std::lock_guard<SpinLock> lock(shard.mtx);
		const size_t n = shard.size.load(std::memory_order_relaxed);
		if (n == shard.buf.size())
			shard.grow();
		shard.buf[(shard.head + n) & (shard.buf.size() - 1)] = std::move(extension);
		shard.size.store(n + 1, std::memory_order_relaxed);
	}
	// Moves the front extension of a shard to out. `shard` is the one tried first, the others follow in cyclic order.
	// It is set to the shard popped from, which the caller should pass to the next call. Returns false if no extension
	// is available at the moment: empty shards are detected from their size without taking their lock, and shards locked
	// by another thread are skipped rather than waited for. So polling an empty queue takes no lock, and the threads
	// popping from a queue holding few extensions do not pile up at the lock of the same shard. A false result only
	// means that the queue is empty if no other thread is using it.
	bool pop(int& shard, Extension& out) {
		return pop(shard, &out, 1) == 1;
	}
	// Moves up to n extensions from the front of the shards to out and returns their number. `shard` is used like in
	// the function above. Less than n extensions are returned if the others are not available at the moment (see above).
	size_t pop(int& shard, Extension* out, size_t n) {
		size_t popped = 0;
		int s = shard % shard_count_;
		for (int k = 0; k < shard_count_ && popped < n; ++k, s = s + 1 < shard_count_ ? s + 1 : 0) {
			Shard& sh = shards_[s];
			if (sh.size.load(std::memory_order_relaxed) == 0)
				continue;
			std::unique_lock<SpinLock> lock(sh.mtx, std::try_to_lock);
			if (!lock.owns_lock())
				continue;
			size_t size = sh.size.load(std::memory_order_relaxed);
			if (size == 0)
				continue;
			for (; size > 0 && popped < n; --size, ++popped) {
				out[popped] = std::move(sh.buf[sh.head]);
				sh.head = (sh.head + 1) & (sh.buf.size() - 1);
			}
			sh.size.store(size, std::memory_order_relaxed);
			shard = s;
		}
		return popped;
	}
	// Number of extensions held. Read without the locks of the shards, so it is only exact if no thread is pushing
	// or popping.
	size_t size() const {
		size_t n = 0;
		for (int s = 0; s < shard_count_; ++s)
			n += shards_[s].size.load(std::memory_order_relaxed);
		return n;
	}
	const Loc band;                 // band limit of the anchored swipe, at least the band of every extension of the queue

private:

	// Aligned to a cache line, so that threads working on different shards do not share one.
	struct alignas(64) Shard {

		static constexpr size_t INITIAL_CAPACITY = 256;

		// Doubles the capacity of the full buffer, moving the extensions to the front of the new one.
		void grow() {
			const size_t n = buf.size();
			std::vector<Extension> new_buf(n == 0 ? INITIAL_CAPACITY : 2 * n);
			for (size_t k = 0; k < n; ++k)
				new_buf[k] = std::move(buf[(head + k) & (n - 1)]);
			buf.swap(new_buf);
			head = 0;
		}

		SpinLock mtx;                   // guards the members below
		std::atomic<size_t> size{ 0 };  // number of extensions held, written under the lock but readable without it
		size_t head = 0;                // position of the front extension in buf
		std::vector<Extension> buf;     // capacity is 0 or a power of 2
	};

	const int shard_count_;
	const std::unique_ptr<Shard[]> shards_;
};

// Receives the extensions popped from the queue once they have been computed.
using ExtensionCallback = std::function<void(std::vector<Extension>&)>;

// Histogram of the bands of the extensions computed over the whole run of the program, i.e. of the widths of their
// diagonal ranges (see Extension::diag_spread()). Only collected if --band-hist is set.
// Adds a partial histogram, whose element b counts the extensions of band b, to it. Thread-safe.
void add_band_histogram(const std::vector<int64_t>& counts);
// Writes the histogram to the file as TSV lines of band and count, for the bands of nonzero count in ascending order.
void write_band_histogram(const std::string& file);

}

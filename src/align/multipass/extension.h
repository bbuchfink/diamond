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
#include <array>
#include <atomic>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <new>
#include <type_traits>
#include <utility>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
#include "basic/sequence.h"
#include "basic/value.h"
#include "util/hsp/approx_hsp.h"
#include "util/geo/geo.h"
#include "util/memory/mem_profile.h"

namespace MultipassPipeline {

// The two directions an anchor is extended into. Used to index the result slots of TargetPair.
enum Direction { RIGHT = 0, LEFT = 1, DIRECTIONS = 2 };

// A query/target pair whose anchor is extended. The two directions are independent extensions, each of
// which refers to the pair and stores its result in it. The extension computed last merges the two results
// into one alignment.
struct TargetPair {
	TargetPair(BlockId query_id, BlockId target_id, Sequence query, Sequence target, const Anchor& anchor, const int8_t* score_table, bool score_table_compact) :
		query_id(query_id),
		target_id(target_id),
		query(query),
		target(target),
		anchor(anchor),
		score_table(score_table),
		score_table_compact(score_table_compact),
		pending(0),
		score{ 0, 0 },
		query_len{ 0, 0 },
		target_len{ 0, 0 }
	{}
	// Stores the result of the extension into the given direction and returns true if this was the last
	// extension of the pair, i.e. the caller has to merge and output it. The two extensions of a pair may be
	// computed by different threads.
	bool finish(Direction direction, Score score, Loc query_len, Loc target_len) {
		this->score[direction] = score;
		this->query_len[direction] = query_len;
		this->target_len[direction] = target_len;
		return pending.fetch_sub(1, std::memory_order_acq_rel) == 1;
	}
	const BlockId query_id, target_id;
	const Sequence query, target;
	const Anchor anchor;            // its score is the one under the score table of the pair
	// Compositionally adjusted matrix of the pair, held by the buffer holding the pair. nullptr if the matrix is not
	// adjusted. Either the data of a CompactScoreTable (score_table_compact) or a table in the layout of
	// DP::AnchoredSwipe::Options::score_table.
	const int8_t* const score_table;
	const bool score_table_compact;
	std::atomic<int> pending;       // extensions buffered and not yet computed
	// Results of the extensions, indexed by direction. A direction that is not extended stays zero.
	Score score[DIRECTIONS];
	Loc query_len[DIRECTIONS], target_len[DIRECTIONS];
};

// The part of the DP matrix of a query/target pair that the anchored swipe computes for an extension: the sub-matrix
// the extension runs into, starting at its first cell, and the diagonals of it to compute. For LEFT extensions, it
// refers to the reversed sequences, so that every extension runs towards increasing coordinates.
struct ExtensionMatrix {
	// [d_first, d_last] is the diagonal range of the extension in the coordinates of the sub-matrix, which is clipped
	// to the diagonals of the sub-matrix.
	ExtensionMatrix(Loc query_start, Loc query_len, Loc target_len, Loc d_first, Loc d_last) :
		query_start(query_start),
		query_len(query_len),
		d_begin(Geo::clip_diag(d_first, query_len, target_len)),
		d_end(Geo::clip_diag(d_last, query_len, target_len) + 1),
		target_len(std::min(target_len, Geo::j(query_len - 1, d_begin) + 1))
	{}
	Loc band() const {
		return d_end - d_begin;
	}
	Loc query_start;                // first query letter of the sub-matrix
	Loc query_len;                  // query letters of the sub-matrix
	Loc d_begin, d_end;             // diagonals to compute, relative to the sub-matrix
	Loc target_len;                 // target letters of the sub-matrix, without the ones the diagonals do not reach
};

// Extension of an anchor into one direction of a query/target pair as an independent local alignment.
struct Extension {
	Extension(TargetPair* pair, Direction direction, Loc i, Loc j, Loc d_min, Loc d_max) :
		pair(pair),
		i(i),
		j(j),
		d_min(d_min),
		d_max(d_max),
		direction(direction),
		band(matrix().band())
	{}
	ExtensionMatrix matrix() const {
		const Loc qlen = pair->query.length(), tlen = pair->target.length();
		if (direction == RIGHT)
			return ExtensionMatrix(i, qlen - i, tlen - j, Geo::diag_sub_matrix(d_min, i, j), Geo::diag_sub_matrix(d_max - 1, i, j));
		// The first cell and the diagonal range in the matrix of the reversed sequences.
		const Loc ir = qlen - 1 - i, jr = tlen - 1 - j;
		return ExtensionMatrix(ir, i + 1, j + 1, Geo::diag_sub_matrix(Geo::rev_diag(d_max, qlen, tlen), ir, jr),
			Geo::diag_sub_matrix(Geo::rev_diag(d_min, qlen, tlen) - 1, ir, jr));
	}
	TargetPair* pair;               // holds the sequences
	Loc i, j;                       // first cell: end of the anchor (RIGHT), last letter before it (LEFT)
	Loc d_min, d_max;               // diagonal range to compute
	Direction direction;
	// Band the anchored swipe computes the extension in (see ExtensionMatrix::band): the number of diagonals of its
	// range that lie in the sub-matrix it runs into.
	Loc band;
};

/* The scores of a compositionally adjusted matrix that depend on the pair: the ones of the true amino acids and the
   mask letter against each other. All other scores of the 32x32 score table of the pair (see
   DP::AnchoredSwipe::Options::score_table) are the same for every pair, so that the table is stored as these scores
   only and expanded into a table holding the common scores when the pair is computed. */
struct CompactScoreTable {

	// Letters with adjusted scores.
	static constexpr int N = TRUE_AA + 1;
	static_assert(MASK_LETTER >= TRUE_AA && MASK_LETTER < 32, "The mask letter is stored as the last letter.");

	// The row of a 32x32 score table that row r is stored from.
	static int row(int r) {
		return (r < TRUE_AA ? r : (int)MASK_LETTER) << 5;
	}

	void compress(const int8_t* table) {
		for (int r = 0; r < N; ++r) {
			const int8_t* const src = table + row(r);
			int8_t* const dst = data + r * N;
			std::copy(src, src + TRUE_AA, dst);
			dst[TRUE_AA] = src[(int)MASK_LETTER];
		}
	}

	// Writes the scores to a 32x32 score table, whose other scores are kept.
	void expand(int8_t* table) const {
		for (int r = 0; r < N; ++r) {
			const int8_t* const src = data + r * N;
			int8_t* const dst = table + row(r);
			std::copy(src, src + TRUE_AA, dst);
			dst[(int)MASK_LETTER] = src[TRUE_AA];
		}
	}

	int8_t data[N * N];
};

// Append-only storage of objects, which are constructed in place in chunks of a fixed size allocated on demand.
// A chunk is never reallocated, so appending neither copies the objects stored nor invalidates pointers to them,
// and it takes one heap allocation per chunk. The memory is freed by the destructor. Not thread-safe.
template<typename T>
struct ChunkedBuffer {

	static_assert(std::is_trivially_destructible<T>::value, "The objects are not destroyed.");
	static_assert(alignof(T) <= alignof(std::max_align_t), "Over-aligned types are not supported.");

	// Size of a chunk in bytes, the size of a huge page.
	static constexpr size_t CHUNK_BYTES = (size_t)1 << 21;
	// Objects per chunk.
	static constexpr size_t CHUNK_SIZE = CHUNK_BYTES / sizeof(T);
	static_assert(sizeof(T) <= CHUNK_BYTES, "The objects have to fit into a chunk.");

	ChunkedBuffer() :
		last_size_(CHUNK_SIZE)
	{}

	ChunkedBuffer(const ChunkedBuffer&) = delete;
	ChunkedBuffer& operator=(const ChunkedBuffer&) = delete;

	~ChunkedBuffer() {
		for (T* chunk : chunks_)
			free_chunk(chunk);
	}

	template<typename... Args>
	T& emplace_back(Args&&... args) {
		if (last_size_ == CHUNK_SIZE) {
			chunks_.push_back(nullptr);
			// A buffer that needs more than one chunk is likely to grow large.
			chunks_.back() = allocate_chunk(chunks_.size() > 1);
			last_size_ = 0;
		}
		T* const p = chunks_.back() + last_size_;
		new (p) T(std::forward<Args>(args)...);
		++last_size_;
		return *p;
	}

	size_t size() const {
		return chunks_.empty() ? 0 : (chunks_.size() - 1) * CHUNK_SIZE + last_size_;
	}

	// Calls f(T* begin, T* end) for the objects of every chunk.
	template<typename F>
	void for_each_chunk(F f) const {
		for (size_t i = 0; i < chunks_.size(); ++i)
			f(chunks_[i], chunks_[i] + (i + 1 < chunks_.size() ? CHUNK_SIZE : last_size_));
	}

private:

#ifdef __linux__
	/* The chunks are mapped directly, aligned to the size of a huge page, and the system is asked to back them by
	   huge pages if huge_pages is set. The buffers of a large run take most of its memory, and with regular pages,
	   mapping them page by page and returning them to the system (which the threads of a process cannot do in
	   parallel) takes a considerable part of its run time. */
	static T* allocate_chunk(bool huge_pages) {
		void* const p = mmap(nullptr, 2 * CHUNK_BYTES, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
		if (p == MAP_FAILED)
			throw std::bad_alloc();
		// Unmap the parts before and after the aligned chunk.
		const uintptr_t begin = (uintptr_t)p, end = begin + 2 * CHUNK_BYTES, chunk = (begin + CHUNK_BYTES - 1) & ~(uintptr_t)(CHUNK_BYTES - 1);
		if (chunk > begin)
			munmap((void*)begin, chunk - begin);
		if (end > chunk + CHUNK_BYTES)
			munmap((void*)(chunk + CHUNK_BYTES), end - (chunk + CHUNK_BYTES));
		if (huge_pages)
			madvise((void*)chunk, CHUNK_BYTES, MADV_HUGEPAGE);
		MEM_ADD("align/multipass-buffers", CHUNK_BYTES);
		return reinterpret_cast<T*>(chunk);
	}

	static void free_chunk(T* chunk) {
		munmap(chunk, CHUNK_BYTES);
		MEM_SUB("align/multipass-buffers", CHUNK_BYTES);
	}
#else
	static T* allocate_chunk(bool) {
		return static_cast<T*>(::operator new(CHUNK_SIZE * sizeof(T)));
	}

	static void free_chunk(T* chunk) {
		::operator delete(chunk);
	}
#endif

	std::vector<T*> chunks_;
	size_t last_size_;              // objects in the last chunk, CHUNK_SIZE if there is none

};

// The extensions created by the tasks run by one worker along with the pairs and score tables they refer to,
// kept for the whole run. Only accessed by that worker, so that no lock is needed. Aligned to a cache line,
// so that the workers do not share one.
struct alignas(64) ExtensionBuffer {

	// Score table of a pair that is not stored as a CompactScoreTable (see TargetPair::score_table).
	using ScoreTable = std::array<int8_t, 32 * 32>;

	template<typename... Args>
	void add_extension(Args&&... args) {
		const Extension& e = extensions.emplace_back(std::forward<Args>(args)...);
		// The anchored swipe requires a band of at least one diagonal.
		assert(e.band >= 1);
		const size_t band = (size_t)e.band;
		if (band >= band_counts.size())
			band_counts.resize(band + 1, 0);
		++band_counts[band];
	}

	/* Stores the score table of a pair (32x32) and returns the pointer for TargetPair::score_table, and whether it
	   was stored as a CompactScoreTable. That is the case if the table equals common_scores in all scores that a
	   CompactScoreTable does not hold. */
	std::pair<const int8_t*, bool> add_score_table(const int8_t* table, const int8_t* common_scores) {
		CompactScoreTable compact;
		compact.compress(table);
		ScoreTable expanded;
		std::copy(common_scores, common_scores + expanded.size(), expanded.begin());
		compact.expand(expanded.data());
		if (std::equal(expanded.begin(), expanded.end(), table))
			return { compact_score_tables.emplace_back(compact).data, true };
		ScoreTable& full = score_tables.emplace_back();
		std::copy(table, table + full.size(), full.begin());
		return { full.data(), false };
	}

	ChunkedBuffer<Extension> extensions;   // appended to by add_extension
	ChunkedBuffer<TargetPair> pairs;
	ChunkedBuffer<CompactScoreTable> compact_score_tables;
	ChunkedBuffer<ScoreTable> score_tables;
	std::vector<size_t> band_counts;       // element b is the number of extensions of band b
};

}

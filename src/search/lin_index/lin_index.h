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
#include <stdint.h>
#include <algorithm>
#include <atomic>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>
#include "basic/value.h"
#include "util/intrin.h"

struct Block;
enum class SeedEncoding;

namespace Search {

struct Config;

/* Seed index of a single sequence block, used by the linear clustering rounds with
   unidirectional coverage. For every seed it holds the single occurrence in the
   longest sequence of the block containing that seed (the pivot), which is the only
   entry the linear stage 1 kernels look at.

   The hash table itself is never stored. What is written to disk once per block
   is only the sorted list of pivot positions, delta encoded as varints, which is
   an order of magnitude smaller than the table. The table is rebuilt from it in
   memory whenever the block takes part in a block combination as the
   representative (reference) block: the seed of a pivot is recomputed from the
   sequence data at its position, so no keys need to be stored. The positions are
   grouped into chunks of a fixed number of entries, each of which starts with an
   absolute position, so that the file can be decoded by all threads in parallel.

   Each seed shape has its own independent section in the file and its own hash
   table. Only one of them is ever held in memory: the shapes are built, and later
   rebuilt and scanned, strictly one after the other, so that the peak memory of
   the index does not grow with the number of shapes. */

/* Open addressing hash table of the index, parameterized by the slot layout. Two
   layouts exist, see PosSlots and LenSlots: the initial build has to pick the pivot
   out of all occurrences of a seed and therefore needs the sequence lengths, while
   the table that is rebuilt from the file only ever holds pivots and gets by with a
   position alone. Both are lock free and their outcome does not depend on the order
   in which the threads insert.

   Common to both is that the position of an entry occupies the low POS_BITS bits of a
   64 bit word and that a seed is identified by a fingerprint of its hash value taken
   from the low bits, which the bucket index does not use. */
template<typename Slots>
struct LinTable {

	static constexpr int POS_BITS = Slots::POS_BITS;
	static constexpr uint64_t POS_MASK = Slots::POS_MASK;
	// Bytes of memory per slot.
	static constexpr int64_t SLOT_BYTES = Slots::SLOT_BYTES;

	void init(const uint64_t size, const int threads) {
		free();
		size_ = size;
		slots_.alloc(size, std::max(threads, 1));
	}

	void free() {
		slots_.free();
		size_ = 0;
	}

	bool empty() const {
		return size_ == 0;
	}

	uint64_t size() const {
		return size_;
	}

	int64_t mem_size() const {
		return (int64_t)size_ * SLOT_BYTES;
	}

	/* The raw slot words, for extracting the entries of the finished table in place.
	   Empty slots are zero and an entry holds its position in the low POS_BITS bits, so
	   sorting this array by those bits collects the entries into one ascending run,
	   which saves materializing them somewhere else. Sorting it invalidates the table:
	   the words no longer belong to the slots they are stored in, so nothing but free()
	   may be called afterwards. */
	uint64_t* words() {
		return slots_.words();
	}

	/* Inserts an occurrence of a seed, keeping the best of all occurrences that map to
	   the same slot, which is defined by the slot layout. len is the length of the
	   sequence the seed comes from and is ignored by layouts that do not store it.
	   Returns true if a new slot was claimed. */
	bool insert(const uint64_t key, const uint64_t pos, const Loc len) {
		const uint64_t h = hash(key), f = Slots::fingerprint(h), value = Slots::value(f, pos, len);
		uint64_t i = bucket(h, size_), probes = 0;
		for (;;) {
			switch (slots_.insert(i, f, value)) {
			case Slots::Result::CLAIMED:
				return true;
			case Slots::Result::UPDATED:
				return false;
			default:
				if (++i == size_)
					i = 0;
				if (++probes == size_)
					throw std::runtime_error("Seed index hash table overflow.");
			}
		}
	}

	// Position of the entry of the seed in the raw sequence data of the indexed block,
	// 0 if the seed is not contained in the table.
	uint64_t operator()(const uint64_t key) const {
		const uint64_t h = hash(key), f = Slots::fingerprint(h);
		uint64_t i = bucket(h, size_), pos;
		for (;;) {
			switch (slots_.probe(i, f, pos)) {
			case Slots::Result::EMPTY:
				return 0;
			case Slots::Result::MATCH:
				return pos;
			default:
				if (++i == size_)
					i = 0;
			}
		}
	}

	static uint64_t hash(uint64_t x) {
		// splitmix64 finalizer
		x += 0x9e3779b97f4a7c15ull;
		x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
		x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
		return x ^ (x >> 31);
	}

	// Maps a hash value to [0, n) without a division (Lemire's fastrange). Only the
	// high bits of the hash value are used, so the fingerprints can be taken from the
	// low ones.
	static uint64_t bucket(const uint64_t h, const uint64_t n) {
		return mul_hi(h, n);
	}

private:

	Slots slots_;
	uint64_t size_ = 0;

};

// Zeroes an array of atomics in parallel.
template<typename T>
static void clear_atomics(std::atomic<T>* v, const uint64_t n, const int threads) {
	std::vector<std::thread> workers;
	for (int i = 0; i < threads; ++i)
		workers.emplace_back([v, n, threads, i] {
			const uint64_t begin = n / threads * i, end = i == threads - 1 ? n : n / threads * (i + 1);
			for (uint64_t j = begin; j < end; ++j)
				v[j].store(T(), std::memory_order_relaxed);
			});
	for (auto& t : workers)
		t.join();
}

// Constants and the outcome of a slot operation, common to all slot layouts.
struct LinSlots {
	enum class Result { EMPTY, MATCH, MISMATCH, CLAIMED, UPDATED };
	// A position in the raw sequence data of the block, occupying the low bits of the
	// slot word. Limits the block to one terabyte of sequence data.
	static constexpr int POS_BITS = 40;
	static constexpr uint64_t POS_MASK = (uint64_t(1) << POS_BITS) - 1;
};

/* Slot layout of the table that is rebuilt from the index file: a single 64 bit word
   holding a 24 bit fingerprint of the seed and the 40 bit position of the entry. An
   empty slot is zero, which is never a valid entry since position 0 lies in the
   leading padding of a sequence set. The file holds nothing but pivots, so no
   occurrence of a seed has to be selected here; of two distinct seeds that collide on
   their fingerprint the one with the smaller position wins, which is arbitrary but
   deterministic. */
struct PosSlots : public LinSlots {

	static constexpr int64_t SLOT_BYTES = sizeof(uint64_t);
	static constexpr int FP_BITS = 64 - POS_BITS;

	static uint64_t fingerprint(const uint64_t h) {
		return std::max(h & ((uint64_t(1) << FP_BITS) - 1), uint64_t(1));
	}

	static uint64_t value(const uint64_t f, const uint64_t pos, const Loc) {
		return (f << POS_BITS) | pos;
	}

	void alloc(const uint64_t n, const int threads) {
		word_.reset(new std::atomic<uint64_t>[n]);
		clear_atomics(word_.get(), n, threads);
	}

	void free() {
		word_.reset();
	}

	uint64_t* words() {
		static_assert(sizeof(std::atomic<uint64_t>) == sizeof(uint64_t), "");
		return reinterpret_cast<uint64_t*>(word_.get());
	}

	Result probe(const uint64_t i, const uint64_t f, uint64_t& pos) const {
		const uint64_t e = word_[i].load(std::memory_order_relaxed);
		if (e == 0)
			return Result::EMPTY;
		if ((e >> POS_BITS) != f)
			return Result::MISMATCH;
		pos = e & POS_MASK;
		return Result::MATCH;
	}

	Result insert(const uint64_t i, const uint64_t f, const uint64_t value) {
		for (;;) {
			uint64_t cur = word_[i].load(std::memory_order_relaxed);
			if (cur == 0) {
				if (word_[i].compare_exchange_weak(cur, value, std::memory_order_relaxed, std::memory_order_relaxed))
					return Result::CLAIMED;
				continue;
			}
			if ((cur >> POS_BITS) != f)
				return Result::MISMATCH;
			// The fingerprints are equal, so comparing the words compares the positions.
			if (cur <= value)
				return Result::UPDATED;
			if (word_[i].compare_exchange_weak(cur, value, std::memory_order_relaxed, std::memory_order_relaxed))
				return Result::UPDATED;
		}
	}

private:

	std::unique_ptr<std::atomic<uint64_t>[]> word_;

};

/* Slot layout of the initial build, which has to find the occurrence of a seed in the
   longest sequence of the block. The length of the sequence a seed comes from is kept
   in the slot, so that the block does not have to be sorted by decreasing length for
   the pivot to be identifiable by its position alone.

   A slot is a 64 bit word holding the length of the sequence complemented to 16 bits
   and the 40 bit position of the entry, plus a 16 bit fingerprint of the seed in a
   separate array. The pivot is then simply the smallest word of a slot: the
   complemented length orders the longest sequence first, and the position breaks ties
   between sequences of equal length. Lengths above LEN_MAX saturate, so seeds of very
   long sequences are ordered by position among themselves.

   Both fields have to be selected in one atomic step for the outcome to be independent
   of the thread order, which is why the length lives in the word and the fingerprint,
   which is fixed once a slot has been claimed, in the extra two bytes rather than the
   other way around. A 16 bit fingerprint lets distinct seeds share a slot with a
   probability of 2^-16 per probe, which costs a negligible fraction of the pivots. */
struct LenSlots : public LinSlots {

	static constexpr int64_t SLOT_BYTES = sizeof(uint64_t) + sizeof(uint16_t);
	static constexpr int FP_BITS = 16;
	static constexpr Loc LEN_MAX = 0xffff;

	static uint64_t fingerprint(const uint64_t h) {
		return std::max(h & ((uint64_t(1) << FP_BITS) - 1), uint64_t(1));
	}

	static uint64_t value(const uint64_t, const uint64_t pos, const Loc len) {
		return ((uint64_t)(LEN_MAX - std::min(len, LEN_MAX)) << POS_BITS) | pos;
	}

	void alloc(const uint64_t n, const int threads) {
		word_.reset(new std::atomic<uint64_t>[n]);
		fp_.reset(new std::atomic<uint16_t>[n]);
		clear_atomics(word_.get(), n, threads);
		clear_atomics(fp_.get(), n, threads);
	}

	void free() {
		word_.reset();
		fp_.reset();
	}

	uint64_t* words() {
		static_assert(sizeof(std::atomic<uint64_t>) == sizeof(uint64_t), "");
		return reinterpret_cast<uint64_t*>(word_.get());
	}

	Result probe(const uint64_t i, const uint64_t f, uint64_t& pos) const {
		const uint64_t cur = fp_[i].load(std::memory_order_relaxed);
		if (cur == 0)
			return Result::EMPTY;
		if (cur != f)
			return Result::MISMATCH;
		pos = word_[i].load(std::memory_order_relaxed) & POS_MASK;
		return Result::MATCH;
	}

	/* Claims the slot for the fingerprint if it is empty, then reduces the word to the
	   minimum of itself and the new value. A slot that has just been claimed still has
	   a zero word, which stands for an unset value and loses against everything. */
	Result insert(const uint64_t i, const uint64_t f, const uint64_t value) {
		uint16_t cur_f = fp_[i].load(std::memory_order_relaxed);
		bool claimed = false;
		if (cur_f == 0) {
			claimed = fp_[i].compare_exchange_strong(cur_f, (uint16_t)f, std::memory_order_relaxed, std::memory_order_relaxed);
			if (claimed)
				cur_f = (uint16_t)f;
		}
		if (cur_f != f)
			return Result::MISMATCH;
		for (;;) {
			uint64_t cur = word_[i].load(std::memory_order_relaxed);
			if (cur != 0 && cur <= value)
				break;
			if (word_[i].compare_exchange_weak(cur, value, std::memory_order_relaxed, std::memory_order_relaxed))
				break;
		}
		return claimed ? Result::CLAIMED : Result::UPDATED;
	}

private:

	std::unique_ptr<std::atomic<uint64_t>[]> word_;
	std::unique_ptr<std::atomic<uint16_t>[]> fp_;

};

struct LinIndex {

	// The table of the file is the one without the sequence lengths: it only ever holds
	// the pivots that the initial build has already selected.
	using Table = LinTable<PosSlots>;

	static constexpr uint64_t MAGIC = 0x786469646e696c64ull;
	static constexpr uint32_t VERSION = 4;
	static constexpr uint64_t POS_MASK = Table::POS_MASK;
	// Load factor of the hash tables.
	static constexpr double HASH_TABLE_LOAD = 0.8;
	static constexpr uint64_t CHUNK_ENTRIES = 4096;
	/* Safety margin applied to the estimate of the number of distinct seeds that sizes the
	   build table. The estimate comes from HyperLogLog counters of a few hundred to a
	   thousand registers, whose relative error is a few percent, and a table too small to
	   hold all distinct seeds overflows. */
	static constexpr double BUILD_TABLE_MARGIN = 1.07;

	// Number of slots needed to hold n entries at the load factor above.
	static uint64_t table_size(const uint64_t n) {
		return std::max<uint64_t>(16, (uint64_t)((double)n / HASH_TABLE_LOAD) + 1);
	}

	// Number of slots of the table of the initial build, which is sized from an estimate
	// of the number of distinct seeds rather than from an exact count.
	static uint64_t build_table_size(const uint64_t expected_pivots) {
		return table_size((uint64_t)((double)expected_pivots * BUILD_TABLE_MARGIN));
	}

	/* Memory of the seed index of a block, while it is being built and while it is
	   resident for the scan of a member block. The two differ by a factor of 1.5625: the
	   build table keeps the sequence lengths in an extra two bytes per slot and is sized
	   from an estimate plus BUILD_TABLE_MARGIN, while the table rebuilt from the file
	   holds a bare position per slot and its size is known exactly. The memory model of
	   the superblock partitioning (combine_blocks.cpp) has to charge whichever of the two
	   is resident at the point it is modelling, so both live here. */
	static uint64_t build_memory(const uint64_t expected_pivots) {
		return build_table_size(expected_pivots) * (uint64_t)LenSlots::SLOT_BYTES;
	}

	static uint64_t index_memory(const uint64_t pivots) {
		return table_size(pivots) * (uint64_t)Table::SLOT_BYTES;
	}

	struct Header {
		uint64_t magic;
		uint32_t version;
		uint32_t shape_count;
		uint64_t seq_count;
		uint64_t raw_len;
		// The seed encoding the positions were enumerated with, see encoding().
		uint64_t seed_encoding;
	};

	/* Seed encoding of the index. The hashed encoding is the preferred one: it computes
	   the seed of a window with a mask and a hash instead of gathering the letters of the
	   shape one by one, which makes the cost per seed independent of the weight of the
	   shape. Its hash collisions are harmless here, since a colliding seed only costs a
	   pivot lookup whose hit the ungapped filters of stage 2 discard again. It cannot hold
	   every shape and does not support the minimizer though, so the letterwise encoding
	   stays as the fallback. Both sides of a block combination have to agree on it, which
	   the header of the file is checked for. */
	static SeedEncoding encoding(const Config& cfg);

	// Directory entry of the section of one seed shape. The section consists of the
	// chunk offset table followed by the varint encoded positions and starts at
	// file_offset.
	struct ShapeIndex {
		uint64_t table_size;
		uint64_t entry_count;
		uint64_t data_size;
		uint64_t file_offset;
	};

	// Reads the header and the shape directory of the file. No hash table is built
	// yet, use build_shape for that. The block has to be hard masked in the same way
	// as when the file was written and has to stay alive for the lifetime of this
	// object.
	LinIndex(const std::string& file_name, const Block& block, const Config& cfg, int threads);
	~LinIndex();

	// Rebuilds the hash table of one seed shape, discarding the table of the shape
	// that was built before.
	void build_shape(int shape_id);
	void free_shape();

	// Position of the pivot occurrence of the seed in the raw sequence data of
	// the indexed block, 0 if the seed is not contained in the index. Refers to the
	// shape that build_shape was last called with.
	uint64_t operator()(const uint64_t key) const {
		return table_(key);
	}

	const Header& header() const {
		return header_;
	}

	// Seed encoding the member block has to be enumerated with to match this index.
	SeedEncoding seed_encoding() const {
		return (SeedEncoding)header_.seed_encoding;
	}

	/* Verifies that the round is configured to enumerate seeds the way the index was
	   built. Not done when the file is opened, which happens before the search
	   configuration of the round is complete. */
	void check_encoding(const Config& cfg) const;

	// Total number of pivots over all seed shapes.
	uint64_t entry_count() const {
		uint64_t n = 0;
		for (const ShapeIndex& s : shape_index_)
			n += s.entry_count;
		return n;
	}

	uint64_t shape_entry_count(int shape_id) const {
		return shape_index_[shape_id].entry_count;
	}

	// Peak size of the hash table in bytes, i.e. that of the largest shape, since
	// only one shape is held in memory at a time.
	int64_t size() const {
		uint64_t n = 0;
		for (const ShapeIndex& s : shape_index_)
			n = std::max(n, s.table_size);
		return (int64_t)n * Table::SLOT_BYTES;
	}

	int64_t table_size() const {
		return table_.mem_size();
	}

	uint64_t insert_count() const {
		return insert_count_;
	}

private:

	const std::string file_name_;
	const Block& block_;
	const int threads_;
	Header header_;
	std::vector<ShapeIndex> shape_index_;
	Table table_;
	uint64_t insert_count_ = 0;
	int current_shape_ = -1;

};

/* Computes the seed positions of the block and writes them to disk, one seed shape
   after the other. The block has to be hard masked in the same way as at search time.

   expected_pivots is an estimate of the number of distinct seeds of the block, maximized
   over the seed shapes, which is what the table of the initial build has to hold. It is
   required: sizing that table from the estimate is what makes it possible to build the
   index in a single pass over the seeds of the block. The size of the table that is
   rebuilt from the file does not depend on it: the number of pivots is known exactly once
   they have been extracted. */
void build_lin_index(Block& block, const std::string& file_name, const Config& cfg, int threads, uint64_t expected_pivots);

// Scans the member block (cfg.target) against the index of the reference block
// (cfg.query) and writes the resulting seed hits to the hit buffer. Iterates over the
// seed shapes, rebuilding the hash table of the reference block for each of them.
void scan_lin_index(Config& cfg);

}
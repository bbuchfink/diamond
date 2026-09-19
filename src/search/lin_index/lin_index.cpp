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

#define NOMINMAX
#include <algorithm>
#include <atomic>
#include <cstring>
#include <thread>
#include <vector>
#include "lin_index.h"
#define _REENTRANT
#include "lib/ips4o/ips4o.hpp"
#include "basic/reduction.h"
#include "basic/shape_config.h"
#include "data/block/block.h"
#include "run/config.h"
#include "search/seed_array/enum_seeds.h"
#include "search/seed_complexity.h"
#include "util/algo/varint.h"
#include "util/io/file.h"
#include "util/log_stream.h"
#include "util/ptr_vector.h"

using std::atomic;
using std::endl;
using std::runtime_error;
using std::string;
using std::to_string;
using std::unique_ptr;
using std::vector;
using std::pair;

namespace Search {

/* Applies the reduction to the raw sequence data on the fly, so that the seed at a
   position can be computed without materializing the reduced sequence. */
struct ReducedSeq {
	Letter operator[](const int i) const {
		return (Letter)Reduction::get_reduction()(letter_mask(ptr[i]));
	}
	const Letter* ptr;
};

// Field width of the hashed seed encoding, which is only instantiated for the four bit
// reductions (see enum_seeds_worker).
static constexpr uint64_t HASHED_SEED_BITS = 4;

// Recomputes the seed of a stored pivot position, mirroring what the seed iterator of the
// encoding yields for that position in enum_seeds.
static bool seed_at(const Letter* seq, const Shape& shape, const SeedEncoding code, uint64_t& key) {
	if (code == SeedEncoding::HASHED)
		return hashed_seed<HASHED_SEED_BITS>(seq, shape, key);
	return shape.set_seed_reduced(key, ReducedSeq{ seq });
}

SeedEncoding LinIndex::encoding(const Config& cfg) {
	/* The minimizer is only implemented over the letterwise encoding, and a shape whose
	   span does not fit into the rolling window of the hashed encoding cannot be masked
	   out of it. */
	return cfg.minimizer_window == 0 && hashed_seeds_supported() ? SeedEncoding::HASHED : SeedEncoding::SPACED_FACTOR;
}

/* The positions of the index are only meaningful together with the encoding they were
   enumerated with, since the seed of a pivot is recomputed from its position. */
void LinIndex::check_encoding(const Config& cfg) const {
	if (header_.seed_encoding != (uint64_t)encoding(cfg))
		throw runtime_error("Seed index was built with a different seed encoding: " + file_name_);
}

// The table of the initial build, which selects the pivots by sequence length.
using BuildTable = LinTable<LenSlots>;

struct InsertCallback {

	InsertCallback(BuildTable& table, const SequenceSet& seqs, const Shape& shape, double seed_complexity_cut) :
		table(table),
		seqs(seqs),
		shape(shape),
		seed_complexity_cut(seed_complexity_cut)
	{}

	bool operator()(const uint64_t key, const uint64_t pos, const uint32_t block_id, uint64_t) {
		/*if (!seed_is_complex(seqs.data(pos), shape, seed_complexity_cut)) {
			++low_complexity;
			return true;
		}*/
		if (table.insert(key, pos, seqs.length((BlockId)block_id)))
			++inserted;
		return true;
	}

	void finish() {}

	BuildTable& table;
	const SequenceSet& seqs;
	const Shape& shape;
	const double seed_complexity_cut;
	uint64_t inserted = 0, low_complexity = 0;

};

/* Orders the slot words of the finished table by the position of their entry. The
   position occupies the low bits of a word, the length of the sequence the entry comes
   from the high ones, so the words cannot be compared as they are. */
struct CmpPos {
	bool operator()(const uint64_t a, const uint64_t b) const {
		return (a & BuildTable::POS_MASK) < (b & BuildTable::POS_MASK);
	}
};

/* Extracts the pivot positions of the finished table, sorted ascending, by sorting the
   slots in place. Every occupied slot holds exactly one pivot, so no lookup is needed to
   tell the pivots apart from the other occurrences of their seed, and the sort is done on
   the table itself, so the positions do not have to be collected anywhere. Empty slots are
   zero and no entry has position zero, so they all sort ahead of the entries, which end up
   as one contiguous run at the end of the array. */
static pair<const uint64_t*, const uint64_t*> extract_positions(BuildTable& table, const int threads) {
	uint64_t* const begin = table.words(), * const end = begin + table.size();
	ips4o::parallel::sort(begin, end, CmpPos(), threads);
	return { std::upper_bound(begin, end, (uint64_t)0, CmpPos()), end };
}

// The longest varint encoding of a 64 bit value takes 10 bytes.
static const size_t MAX_VARINT_LEN = 10;
// Size of the buffer the delta encoded positions are streamed through.
static const size_t POSITION_BUFFER = 4 * MEGABYTES;

/* Delta encodes the sorted pivot positions into chunks of CHUNK_ENTRIES entries and
   appends them to the file. Each chunk starts with an absolute position and can
   therefore be decoded independently of the others.

   There is one position per distinct seed of the block and the build table they are read
   from is still resident here, so they are streamed through a fixed size buffer instead
   of being collected in memory. That leaves the chunk offset table, which precedes them
   in the section: its size follows from the number of entries and is therefore known up
   front, but its content only becomes known as the positions are encoded. So a
   placeholder is put down and overwritten at the end, the same way build_lin_index
   handles the shape directory. */
static void write_positions(const uint64_t* i, const uint64_t* const end, File& out, LinIndex::ShapeIndex& shape_index) {
	const uint64_t entries = (uint64_t)(end - i);
	vector<uint64_t> chunk_offset((size_t)((entries + LinIndex::CHUNK_ENTRIES - 1) / LinIndex::CHUNK_ENTRIES), 0);
	const size_t offset_bytes = chunk_offset.size() * sizeof(uint64_t);

	shape_index.entry_count = entries;
	shape_index.file_offset = (uint64_t)out.tell();
	if (offset_bytes)
		out.write(chunk_offset.data(), offset_bytes);

	vector<char> buf(POSITION_BUFFER);
	char* const buf_begin = buf.data(), * const flush_at = buf_begin + buf.size() - MAX_VARINT_LEN;
	char* p = buf_begin;
	uint64_t prev = 0, flushed = 0, count = 0;
	for (; i < end; ++i) {
		if (p > flush_at) {
			out.write(buf_begin, (size_t)(p - buf_begin));
			flushed += (uint64_t)(p - buf_begin);
			p = buf_begin;
		}
		// The words of the sorted table still carry the sequence length in their high bits.
		const uint64_t pos = *i & BuildTable::POS_MASK;
		if (count % LinIndex::CHUNK_ENTRIES == 0) {
			chunk_offset[(size_t)(count / LinIndex::CHUNK_ENTRIES)] = flushed + (uint64_t)(p - buf_begin);
			prev = 0;
		}
		else if (pos <= prev)
			throw runtime_error("Seed positions are not sorted.");
		p = write_varuint64(pos - prev, p);
		prev = pos;
		++count;
	}
	out.write(buf_begin, (size_t)(p - buf_begin));
	const uint64_t n = flushed + (uint64_t)(p - buf_begin);
	shape_index.data_size = n;

	if (offset_bytes) {
		const int64_t section_end = out.tell();
		out.seek((int64_t)shape_index.file_offset);
		out.write(chunk_offset.data(), offset_bytes);
		out.seek(section_end);
	}
	*log_stream << "Seed index: section size=" << offset_bytes + n
		<< " bytes per entry=" << (entries ? (double)n / entries : 0.0) << endl;
}

// Builds the hash table of one seed shape, extracts its pivot positions and appends
// them to the file. The table is released at the end, so that only one shape occupies
// memory at a time.
static LinIndex::ShapeIndex build_shape_index(Block& block, const int shape_id, File& out, const Config& cfg, const int threads, const uint64_t expected_pivots) {
	const SequenceSet& seqs = block.seqs();
	const Shape& shape = shapes[shape_id];
	const string label = " (shape " + to_string(shape_id) + ")";

	TaskTimer timer;
	const auto partition = seqs.partition(threads);
	/* filter_masked_seeds: a masked letter on a match position of the shape has to
	   invalidate the seed of its window, which the hashed encoding only does if it is
	   asked for. The letterwise encoding always does it and ignores the flag. */
	const EnumCfg enum_cfg{ &partition, shape_id, shape_id + 1, LinIndex::encoding(cfg), nullptr, true, false, cfg.seed_complexity_cut,
		cfg.soft_masking, cfg.minimizer_window, false, false, cfg.sketch_size, nullptr };

	LinIndex::ShapeIndex shape_index;
	/* The build table holds one entry per distinct seed, which is only known exactly once
	   the seeds have been enumerated. It is sized from the estimate of the caller plus a
	   safety margin, which saves a pass over the seeds of the block just to count them. The
	   estimate is the maximum over the seed shapes, so it is on the safe side for all but
	   the most frequent one. */
	const uint64_t build_table_size = LinIndex::build_table_size(expected_pivots);
	*log_stream << "Seed index: shape=" << shape_id << " expected pivots=" << expected_pivots
		<< " table_size=" << build_table_size << " memory=" << (int64_t)build_table_size * BuildTable::SLOT_BYTES << endl;

	timer.go("Allocating seed index" + label);
	BuildTable table;
	table.init(build_table_size, threads);

	timer.go("Building seed index" + label);
	uint64_t inserted = 0, low_complexity = 0;
	{
		PtrVector<InsertCallback> v;
		for (int i = 0; i < threads; ++i)
			v.push_back(new InsertCallback(table, seqs, shape, cfg.seed_complexity_cut));
		enum_seeds(block, v, &no_filter, enum_cfg);
		for (int i = 0; i < threads; ++i) {
			inserted += v[i].inserted;
			low_complexity += v[i].low_complexity;
		}
	}
	*message_stream << "Block size (sequences): " << block.seqs().mem_size() << ", block size (total): " << block.mem_size() << ", index size: "
		<< table.mem_size() << ", slots used: " << (int64_t)inserted * BuildTable::SLOT_BYTES << ", low complexity seeds: " << low_complexity
		<< ", load: " << (double)inserted / build_table_size << endl;

	timer.go("Extracting seed positions" + label);
	const auto positions = extract_positions(table, threads);

	timer.go("Writing seed positions" + label);
	write_positions(positions.first, positions.second, out, shape_index);
	table.free();
	// The table that is rebuilt from the file holds nothing but the pivots that were just
	// written, so its size is known exactly and independently of any estimate.
	shape_index.table_size = LinIndex::table_size(shape_index.entry_count);
	timer.finish();
	return shape_index;
}

void build_lin_index(Block& block, const string& file_name, const Config& cfg, int threads, const uint64_t expected_pivots) {
	if (expected_pivots == 0)
		throw runtime_error("Building the seed index requires an estimate of the number of distinct seeds of the block.");
	const SequenceSet& seqs = block.seqs();
	if ((uint64_t)seqs.raw_len() > LinIndex::POS_MASK)
		throw runtime_error("Block size exceeds the maximum supported by the seed index.");
	threads = std::max(threads, 1);

	LinIndex::Header header;
	header.magic = LinIndex::MAGIC;
	header.version = LinIndex::VERSION;
	header.shape_count = (uint32_t)shapes.count();
	header.seq_count = (uint64_t)seqs.size();
	header.raw_len = (uint64_t)seqs.raw_len();
	header.seed_encoding = (uint64_t)LinIndex::encoding(cfg);
	vector<LinIndex::ShapeIndex> shape_index(shapes.count());

	// The shape directory is only known once all sections have been written, so a
	// placeholder is put down first and overwritten at the end.
	File out(file_name, "wb");
	out.write(header);
	out.write(shape_index.data(), shape_index.size() * sizeof(LinIndex::ShapeIndex));

	for (int i = 0; i < shapes.count(); ++i)
		shape_index[i] = build_shape_index(block, i, out, cfg, threads, expected_pivots);

	out.seek(0);
	out.write(header);
	out.write(shape_index.data(), shape_index.size() * sizeof(LinIndex::ShapeIndex));
	out.close();
}

LinIndex::LinIndex(const string& file_name, const Block& block, const Config& cfg, int threads) :
	file_name_(file_name),
	block_(block),
	threads_(std::max(threads, 1))
{
	const SequenceSet& seqs = block.seqs();
	File in(file_name, "rb");
	in.read(header_);
	if (header_.magic != MAGIC)
		throw runtime_error("Invalid seed index file: " + file_name);
	if (header_.version != VERSION)
		throw runtime_error("Invalid seed index file version: " + file_name);
	if (header_.shape_count != (uint32_t)shapes.count())
		throw runtime_error("Seed index has a different number of shapes: " + file_name);
	if (header_.seq_count != (uint64_t)seqs.size() || header_.raw_len != (uint64_t)seqs.raw_len())
		throw runtime_error("Seed index does not match the block: " + file_name);
	shape_index_.resize(header_.shape_count);
	in.read(shape_index_.data(), shape_index_.size() * sizeof(ShapeIndex));
	for (const ShapeIndex& s : shape_index_)
		// Guarantees that a lookup of a seed that is not in the index terminates on a
		// blank slot instead of probing the table forever.
		if (s.entry_count >= s.table_size)
			throw runtime_error("Invalid seed index file: " + file_name);
	in.close();
}

void LinIndex::free_shape() {
	table_.free();
	current_shape_ = -1;
}

// The section of the file that holds the pivot positions of one seed shape.
struct PivotPositions {
	PivotPositions(const string& file_name, const LinIndex::ShapeIndex& shape_index) :
		entry_count(shape_index.entry_count),
		chunk_offset((shape_index.entry_count + LinIndex::CHUNK_ENTRIES - 1) / LinIndex::CHUNK_ENTRIES),
		data(shape_index.data_size)
	{
		File in(file_name, "rb");
		in.seek((int64_t)shape_index.file_offset);
		in.read(chunk_offset.data(), chunk_offset.size() * sizeof(uint64_t));
		in.read(data.data(), data.size());
		in.close();
	}
	const uint64_t entry_count;
	vector<uint64_t> chunk_offset;
	vector<char> data;
};

/* Decodes the pivot positions in parallel and recomputes the seed at each of them.
   Every chunk starts with an absolute position, so the chunks can be handed out to the
   threads independently. The seeds are inserted into the table and their number is
   returned. */
static uint64_t scan_pivots(const PivotPositions& positions, const SequenceSet& seqs, const Shape& shape, const uint64_t raw_len,
	const SeedEncoding code, LinIndex::Table& table, const string& file_name, const int threads)
{
	const uint64_t chunk_count = (uint64_t)positions.chunk_offset.size(), entry_count = positions.entry_count;
	atomic<uint64_t> next(0), inserts(0);
	vector<std::thread> workers;
	for (int i = 0; i < threads; ++i)
		workers.emplace_back([&] {
			uint64_t key, ins = 0;
			for (uint64_t c = next++; c < chunk_count; c = next++) {
				const uint64_t begin = c * LinIndex::CHUNK_ENTRIES, n = std::min(LinIndex::CHUNK_ENTRIES, entry_count - begin);
				const char* p = positions.data.data() + positions.chunk_offset[c];
				uint64_t pos = 0;
				for (uint64_t j = 0; j < n; ++j) {
					const auto d = read_varuint64(p);
					pos += d.first;
					p = d.second;
					if (pos >= raw_len)
						throw runtime_error("Invalid seed position in file " + file_name);
					if (!seed_at(seqs.data(pos), shape, code, key))
						continue;
					// The file holds nothing but pivots, so no sequence length is
					// needed to select one.
					table.insert(key, pos, 0);
					++ins;
				}
			}
			inserts.fetch_add(ins, std::memory_order_relaxed);
			});
	for (auto& t : workers)
		t.join();
	return inserts.load(std::memory_order_relaxed);
}

void LinIndex::build_shape(const int shape_id) {
	if (shape_id == current_shape_)
		return;
	const ShapeIndex& shape_index = shape_index_[shape_id];
	const string label = " (shape " + to_string(shape_id) + ")";

	TaskTimer timer;
	// The table of the previous shape is released before the positions of the new one
	// are read, so that only one shape occupies memory at a time.
	free_shape();
	timer.go("Loading seed positions" + label);
	const PivotPositions positions(file_name_, shape_index);

	timer.go("Allocating seed index" + label);
	table_.init(shape_index.table_size, threads_);

	timer.go("Rebuilding seed index" + label);
	insert_count_ = scan_pivots(positions, block_.seqs(), shapes[shape_id], header_.raw_len,
		seed_encoding(), table_, file_name_, threads_);
	current_shape_ = shape_id;
	timer.finish();
}

LinIndex::~LinIndex() {
}

}
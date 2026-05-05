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

#include <unordered_map>
#include "multinode.h"
#include <inttypes.h>
#include "volume.h"
#include "data/sequence_file.h"
#include "util/parallel/simple_thread_pool.h"
#include "util/data_structures/queue.h"

using std::ostringstream;
using std::thread;
using std::vector;
using std::string;
using std::unique_ptr;
using std::atomic;
using std::ofstream;
using std::ifstream;
using std::runtime_error;
using std::unordered_set;
using std::endl;
using std::pair;
using std::tuple;

template<typename T>
void exchange_if_smaller(std::atomic<T>& value, T replacement) {
	T current = value.load(std::memory_order_relaxed);
	while (current > replacement &&
		!value.compare_exchange_weak(
			current,
			replacement,
			std::memory_order_relaxed,
			std::memory_order_relaxed)) {}
}

template<typename T>
void exchange_if_larger(std::atomic<T>& value, T replacement) {
	T current = value.load(std::memory_order_relaxed);
	while (current < replacement &&
		!value.compare_exchange_weak(
			current,
			replacement,
			std::memory_order_relaxed,
			std::memory_order_relaxed)) {
	}
}

struct RepChunk {
	static constexpr size_t POISON = std::numeric_limits<size_t>::max();
	size_t output = POISON;
	string data;
	vector<size_t> record_offsets;
	bool operator==(const RepChunk& rhs) const {
		return output == rhs.output;
	}
};

static void flush_buffer(TextBuffer& buffer, vector<size_t>& record_offsets, size_t output, std::unique_ptr<Queue<RepChunk>>& queue) {
	if (buffer.size() == 0)
		return;
	RepChunk chunk;
	chunk.output = output;
	chunk.data.assign(buffer.data(), buffer.size());
	chunk.record_offsets = std::move(record_offsets);
	queue->enqueue(std::move(chunk));
	buffer.clear();
	record_offsets.clear();
}

static tuple<OId, uint64_t, uint64_t> write_reps(Job& job, const VolumedFile& volumes, ptrdiff_t idx, const string& out_dir, bool single_out_file, FileStack& reps_list,
	atomic<OId>& count_all, atomic<OId>& min_all, atomic<OId>& max_all) {
	job.log("Writing representatives. Volume=%lli/%lli Records=%s", idx + 1, volumes.size(), Util::String::format(volumes[idx].record_count).c_str());
	string id_file = job.base_dir() + "rep_ids" + std::to_string(idx);
	ifstream rep_ids(id_file);
	if (!rep_ids)
		throw runtime_error("Error opening file " + id_file);
	string id;
	unordered_set<OId> rep_id_set;
	while (rep_ids >> id) {
		OId rep_id = std::atoll(id.c_str());
		rep_id_set.insert(rep_id);
	}
	rep_ids.close();
	remove_tmp_file(id_file);

	const SequenceFile::Flags flags = SequenceFile::Flags::SEQS | SequenceFile::Flags::TITLES | SequenceFile::Flags::NEED_LETTER_COUNT;
	string out_file = single_out_file ? out_dir + "reps_all.faa" : out_dir + std::to_string(idx) + ".faa";
	const string offset_file = out_file + ".faidx";
	ofstream out(out_file, std::ios::out | std::ios::app | std::ios::binary);
	if (!out)
		throw runtime_error("Error opening file " + out_file);
	unique_ptr<ofstream> offsets;
	if (single_out_file) {
		offsets = std::make_unique<ofstream>(offset_file, std::ios::out | std::ios::app | std::ios::binary);
		if (!offsets || !*offsets)
			throw runtime_error("Error opening file " + offset_file);
	}
	else {
		count_all.store(0, std::memory_order_relaxed);
		min_all.store(std::numeric_limits<OId>::max(), std::memory_order_relaxed);
		max_all.store(0, std::memory_order_relaxed);
	}
	out.seekp(0, std::ios::end);
	const std::ofstream::pos_type out_pos = out.tellp();
	if (!out || out_pos == std::ofstream::pos_type(-1))
		throw runtime_error("Error determining file offset: " + out_file);
	uint64_t out_offset = (uint64_t)out_pos;
	const int formatter_count = std::max(1, config.threads_);
	const size_t flush_size = 64 * 1024;
	unique_ptr<Queue<RepChunk>> queue;
	SimpleThreadPool pool;
	std::thread::id writer_thread;
	queue = std::make_unique<Queue<RepChunk>>(std::max<size_t>(16, formatter_count * 4), formatter_count, 1, RepChunk());
	auto writer = [&](const atomic<bool>& stop) {
		RepChunk chunk;
		while (!stop.load(std::memory_order_relaxed) && queue->wait_and_dequeue(chunk)) {
			out.write(chunk.data.data(), chunk.data.size());
			if (single_out_file) {
				for (const size_t record_offset : chunk.record_offsets)
					(*offsets) << out_offset + record_offset << '\n';
				out_offset += chunk.data.size();
			}
		}
		};
	writer_thread = pool.spawn(writer);

	unique_ptr<SequenceFile> file(new FastaFile({ volumes[idx].path}, flags, amino_acid_traits));
	uint64_t bytes = 0, ms = 0;
	TaskTimer timer;
	Block* b = file->load_seqs(INT64_MAX);
	ms += timer.microseconds();
	bytes += b->raw_bytes();
	const size_t seq_count = b->seqs().size();
	atomic<size_t> next(0);
	atomic<uint64_t> bytes_all(0);
	vector<std::thread::id> formatter_threads;
	atomic<uint64_t> letters_all(0);
	atomic<OId> count_this_volume(0);
	auto worker = [&](const atomic<bool>& stop) {
		TextBuffer buffer;
		vector<size_t> record_offsets;
		size_t j;
		OId count = 0, min = std::numeric_limits<OId>::max(), max = 0;
		uint64_t bytes = 0;
		uint64_t letters = 0;
		while (!stop.load(std::memory_order_relaxed) && (j = next.fetch_add(1, std::memory_order_relaxed), j < seq_count)) {
			const OId oid = std::atoll(b->ids()[j]);
			if (rep_id_set.find(oid) == rep_id_set.end())
				continue;
			record_offsets.push_back(buffer.size());
			Util::Seq::format(b->seqs()[j], b->ids()[j], nullptr, buffer, "fasta", amino_acid_traits);
			if (buffer.size() >= flush_size) {
				bytes += buffer.size();
				flush_buffer(buffer, record_offsets, idx, queue);
			}
			++count;
			letters += b->seqs()[j].length();
			if(oid < min)
				min = oid;
			if(oid > max)
				max = oid;
		}
		bytes += buffer.size();
		flush_buffer(buffer, record_offsets, idx, queue);
		count_all.fetch_add(count, std::memory_order_relaxed);
		count_this_volume.fetch_add(count, std::memory_order_relaxed);
		bytes_all.fetch_add(bytes, std::memory_order_relaxed);
		letters_all.fetch_add(letters, std::memory_order_relaxed);
		exchange_if_smaller(min_all, min);
		exchange_if_larger(max_all, max);
		};
	for (int j = 0; j < formatter_count; ++j)
		formatter_threads.push_back(pool.spawn(worker));
	try {
		pool.join(formatter_threads.begin(), formatter_threads.end());
	}
	catch (...) {
		delete b;
		for (int i = 0; i < formatter_count; ++i)
			queue->close();
		pool.join(writer_thread);
		throw;
	}
	delete b;

	for (int i = 0; i < formatter_count; ++i)
		queue->close();
	pool.join(writer_thread);
	if (!out)
		throw runtime_error("Error writing representative block");
	if (single_out_file && (!offsets || !*offsets))
		throw runtime_error("Error writing representative offset file");
	file.reset();
	remove_tmp_file(volumes[idx].path);
	if (!single_out_file || idx == volumes.size() - 1) {
		ostringstream ss;
		ss << out_file << '\t' << count_all.load() << '\t' << min_all.load() << '\t' << max_all.load() + 1 << endl;
		reps_list.push(ss.str());
	}
	return { count_this_volume, letters_all, bytes_all.load(std::memory_order_relaxed)};
}

pair<string, uint64_t> get_reps(Job& job, const VolumedFile& volumes) {
	if (job.last_round()) {		
		return { string(),0 };
	}
	const string base_dir = job.base_dir() + PATH_SEPARATOR + "reps" + PATH_SEPARATOR, qpath = base_dir + "queue";
	mkdir(base_dir);
	unique_ptr<FileStack> reps_list(new FileStack(base_dir + "reps.tsv"));

	Atomic q(qpath);
	atomic<int> volumes_processed(0);
	OId cluster_count(0);
	uint64_t bytes = 0, letters = 0;
	int64_t v = 0;
	TaskTimer timer;
	const bool single_out_file = !job.last_round() && !ends_with(job.steps().at(job.round() + 1), "_lin");
	atomic<OId> count_all(0);
	atomic<OId> min_all(std::numeric_limits<OId>::max());
	atomic<OId> max_all(0);
	while (v = q.fetch_add(), v < (int64_t)volumes.size()) {
		volumes_processed.fetch_add(1, std::memory_order_relaxed);
		auto [count, seq_letters, bytes_written] = write_reps(job, volumes, v, base_dir, single_out_file, *reps_list, count_all, min_all, max_all);
		cluster_count += count;
		letters += seq_letters;
		bytes += bytes_written;
	}
	q.remove();
	if (!config.fasta_index_file.empty())
		remove_tmp_file(config.fasta_index_file);
	volumes.remove(job.round() > 0, false);	
	job.log("Representatives written: %" PRIu64 " letters: %" PRIu64, cluster_count, letters);
	const int64_t t = timer.microseconds();
	job.log("Wrote %zu bytes to disk at %.2f MB/s", bytes, (double)bytes / MEGABYTES / (t / 1e6));
	Atomic finished(base_dir + "finished");
	finished.fetch_add(volumes_processed);
	finished.await((int)volumes.size());
	finished.remove();
	const string out = reps_list->file_name();
	reps_list.reset();
	return { out, letters };
}

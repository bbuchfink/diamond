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
#include <inttypes.h>
#include "multinode.h"
#include "volume.h"
#include "data/sequence_file.h"
#include "util/parallel/simple_thread_pool.h"
#include "util/data_structures/queue.h"
#include "util/log_stream.h"
#include "util/io/file.h"

using std::ostringstream;
using std::thread;
using std::vector;
using std::string;
using std::unique_ptr;
using std::atomic;
using std::ofstream;
using std::ifstream;
using std::runtime_error;
using std::unordered_map;
using std::endl;
using std::pair;
using std::tuple;

struct RepChunk {
	static constexpr size_t POISON = std::numeric_limits<size_t>::max();
	size_t order = POISON;
	string data;
	vector<size_t> record_offsets;
	bool operator==(const RepChunk& rhs) const {
		return order == rhs.order;
	}
};

struct RepWriteConfig {
	string output_path;
	bool single_out_file;
	FileStack& reps_list;
};

// Wakes up all threads blocked on the queue. The stop flag of the thread pool is only polled in
// between queue operations, so a thread that is waiting inside the queue semaphore will never
// observe it, and joining that thread (e.g. in ~SimpleThreadPool while an exception is unwinding
// the stack) deadlocks.
struct QueueAbortGuard {
	Queue<RepChunk>* queue;
	~QueueAbortGuard() {
		if (queue)
			queue->abort();
	}
};

static tuple<OId, uint64_t, uint64_t> write_reps(Job& job, const VolumedFile& volumes, size_t idx, const RepWriteConfig& cfg,
	const vector<OId>& clustering,
	atomic<OId>& count_all, atomic<OId>& min_all, atomic<OId>& max_all) {
	job.log("Writing representatives. Volume=%lli/%lli Records=%s", idx + 1, volumes.size(), Util::String::format(volumes[idx].record_count).c_str());

	const SequenceFile::Flags flags = SequenceFile::Flags::SEQS | SequenceFile::Flags::TITLES | SequenceFile::Flags::NEED_LETTER_COUNT;
	string out_file = cfg.single_out_file ? cfg.output_path : cfg.output_path + std::to_string(idx) + ".faa";
	const string offset_file = out_file + ".faidx";
	ofstream out(out_file, std::ios::out | (cfg.single_out_file && idx == 0 ? std::ios::trunc : std::ios::app ) | std::ios::binary);
	if (!out)
		throw runtime_error("Error opening file " + out_file);
	unique_ptr<ofstream> offsets;
	if (cfg.single_out_file) {
		offsets.reset(new ofstream(offset_file, std::ios::out | std::ios::app | std::ios::binary));
		if (!offsets || !*offsets)
			throw runtime_error("Error opening file " + offset_file);
	}
	if (!cfg.single_out_file) {
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
	unique_ptr<Queue<RepChunk>> queue;
	SimpleThreadPool pool;
	std::thread::id writer_thread;
	queue.reset(new Queue<RepChunk>(std::max<size_t>(16, formatter_count * 4), formatter_count, 1, RepChunk()));
	// Declared after `pool` and therefore destructed before it: any exception thrown below releases
	// the worker threads from the queue before the pool destructor joins them.
	QueueAbortGuard abort_guard{ queue.get() };
	auto writer = [&](const atomic<bool>& stop) {
		// Releases the formatter threads that are blocked in enqueue() if this thread exits early
		// (stop flag set or exception thrown).
		QueueAbortGuard writer_guard{ queue.get() };
		unordered_map<size_t, RepChunk> pending;
		size_t next_expected = 0;
		RepChunk chunk;
		while (!stop.load(std::memory_order_relaxed) && queue->wait_and_dequeue(chunk)) {
			pending.emplace(chunk.order, std::move(chunk));
			auto it = pending.find(next_expected);
			while (it != pending.end()) {
				out.write(it->second.data.data(), it->second.data.size());
				if (cfg.single_out_file) {
					for (const size_t record_offset : it->second.record_offsets)
						(*offsets) << out_offset + record_offset << '\n';
					out_offset += it->second.data.size();
				}
				pending.erase(it);
				++next_expected;
				it = pending.find(next_expected);
			}
		}
		};
	writer_thread = pool.spawn(writer);

	uint64_t bytes = 0, ms = 0;
	TaskTimer timer;
	unique_ptr<SequenceFile> file;
	Block* b = nullptr;
	try {
		file.reset(new FastaFile({ volumes[idx].path }, flags, amino_acid_traits));
		b = file->load_seqs(INT64_MAX);
	}
	catch (const EmptyFileError& e) {
		remove_tmp_file(volumes[idx].path);
		return std::make_tuple<OId, uint64_t, uint64_t>(0, 0, 0);
	}
	catch (const std::exception& e) {
		throw runtime_error("Error loading representative volume " + volumes[idx].path + "): " + e.what());
	}
	ms += timer.microseconds();
	bytes += b->raw_bytes();
	const size_t seq_count = b->seqs().size();
	const size_t target_chunks = std::max<size_t>(1, (size_t)formatter_count * 16);
	const size_t chunk_size = std::max<size_t>(1, std::min<size_t>(seq_count / target_chunks + 1, (size_t)8192));
	const size_t chunk_count = (seq_count + chunk_size - 1) / chunk_size;
	atomic<size_t> next(0);
	atomic<uint64_t> bytes_all(0);
	vector<std::thread::id> formatter_threads;
	atomic<uint64_t> letters_all(0);
	atomic<OId> count_this_volume(0);
	auto worker = [&](const atomic<bool>& stop) {
		size_t c;
		OId count = 0, min = std::numeric_limits<OId>::max(), max = 0;
		uint64_t bytes = 0;
		uint64_t letters = 0;
		while (!stop.load(std::memory_order_relaxed) && (c = next.fetch_add(1, std::memory_order_relaxed), c < chunk_count)) {
			TextBuffer buffer;
			vector<size_t> record_offsets;
			const size_t j_begin = c * chunk_size;
			const size_t j_end = std::min(j_begin + chunk_size, seq_count);
			for (size_t j = j_begin; j < j_end; ++j) {
				const OId oid = std::atoll(b->ids()[j]);
				if (oid >= clustering.size() || clustering[oid] != oid)
					continue;
				record_offsets.push_back(buffer.size());
				Util::Seq::format(b->seqs()[j], b->ids()[j], nullptr, buffer, "fasta", amino_acid_traits);
				++count;
				letters += b->seqs()[j].length();
				if(oid < min)
					min = oid;
				if(oid > max)
					max = oid;
			}
			bytes += buffer.size();
			RepChunk chunk;
			chunk.order = c;
			chunk.data.assign(buffer.data(), buffer.size());
			chunk.record_offsets = std::move(record_offsets);
			queue->enqueue(std::move(chunk));
		}
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
		queue->abort();
		pool.join(writer_thread);
		throw;
	}
	delete b;

	for (int i = 0; i < formatter_count; ++i)
		queue->close();
	pool.join(writer_thread);
	if (!out)
		throw runtime_error("Error writing representative block");
	if (cfg.single_out_file && (!offsets || !*offsets))
		throw runtime_error("Error writing representative offset file");
	file.reset();
	remove_tmp_file(volumes[idx].path);
	if (!cfg.single_out_file || idx == volumes.size() - 1) {
		ostringstream ss;
		ss << out_file << '\t' << count_all.load() << '\t' << min_all.load() << '\t' << max_all.load() + 1 << endl;
		cfg.reps_list.push(ss.str());
	}
	return std::make_tuple<OId, uint64_t, uint64_t>(count_this_volume, letters_all, bytes_all.load(std::memory_order_relaxed));
}

pair<string, uint64_t> get_reps(Job& job, const string& round_minichunks, unique_ptr<vector<OId>> clustering) {
	if (job.last_round()) {
		VolumedFile volumes(round_minichunks);
		volumes.remove(false, true, false);
		return { string(), 0 };
	}
	const bool single_out_file = !ends_with(job.steps().at(job.round() + 1), "_lin");
	const string base_dir = job.base_dir() + PATH_SEPARATOR + "rep_minichunks" + PATH_SEPARATOR, qpath = base_dir + "queue";
	const string reps_list_name = base_dir + "reps.tsv";
	job.make_temp_dir(base_dir);
	Atomic get_reps_lock(base_dir + "get_reps_lock", job), get_reps_done(base_dir + "get_reps_done", job);
	Atomic finished(base_dir + "finished", job);
	Atomic letter_count(base_dir + "letter_count", job);
	Atomic q(qpath, job);

	if ((!single_out_file || get_reps_lock.fetch_add() == 0) && get_reps_done.get() == 0) {
		if (!clustering) {
			const string path = job.base_dir() + "clusters.bin";
			File in(path, "rb");
			const size_t count = job.max_oid() + 1;
			if ((size_t)in.size() != count * sizeof(OId))
				throw runtime_error("Invalid binary clustering file size: " + path);
			clustering = std::make_unique<vector<OId>>(count);
			in.read(clustering->data(), clustering->size() * sizeof(OId));
			in.close();
		}
		unique_ptr<FileStack> reps_list(new FileStack(reps_list_name));
		
		OId cluster_count(0);
		uint64_t bytes = 0;
		int64_t v = 0;
		TaskTimer timer;

		const RepWriteConfig cfg{
			single_out_file ? base_dir + "reps_all.faa" : base_dir,
			single_out_file,
			*reps_list
		};
		atomic<OId> count_all(0);
		atomic<OId> min_all(std::numeric_limits<OId>::max());
		atomic<OId> max_all(0);
		VolumedFile volumes(round_minichunks);
		while (v = q.fetch_add(), v < (int64_t)volumes.size()) {
			OId count;
			uint64_t seq_letters, bytes_written;
			std::tie(count, seq_letters, bytes_written) = write_reps(job, volumes, v, cfg, *clustering, count_all, min_all, max_all);
			cluster_count += count;
			bytes += bytes_written;
			letter_count.fetch_add(seq_letters);
			finished.fetch_add();
		}
		finished.await((int)volumes.size());		
		if (!config.fasta_index_file.empty())
			remove_tmp_file(config.fasta_index_file);
		volumes.remove(false, true, false);
		job.log("Representatives written: %" PRIu64 " letters: %" PRIu64, cluster_count, letter_count.get());
		const int64_t t = timer.microseconds();
		job.log("Wrote %zu bytes to disk at %.2f MB/s", bytes, (double)bytes / MEGABYTES / (t / 1e6));
		reps_list.reset();
		get_reps_done.fetch_add();
		job.finish_step();
	}
	else {
		get_reps_done.await(1);
	}
	return { reps_list_name, letter_count.get() };
}

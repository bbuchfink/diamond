/****
DIAMOND protein aligner
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
#include <chrono>
#include <condition_variable>
#include <exception>
#include <mutex>
#include <thread>
#include <unordered_map>
#include "basic/config.h"
#include "util/string/string.h"
#include "util/parallel/filestack.h"
#include "util/parallel/atomic.h"
#include "util/system/system.h"
#include "volume.h"
#include "util/algo/hyperloglog.h"

struct Block;

struct ClusterStats {
	uint64_t hits_evalue_filtered = 0, extensions_computed = 0, hits_filtered = 0, seeds_considered = 0, seeds_indexed = 0;
	MaskingStat masking_stat;
	std::mutex mtx;
	void add(const ClusterStats& s) {
		std::lock_guard<std::mutex> lock(mtx);
		hits_evalue_filtered += s.hits_evalue_filtered;
		extensions_computed += s.extensions_computed;
		hits_filtered += s.hits_filtered;
		masking_stat += s.masking_stat;
		seeds_considered += s.seeds_considered;
		seeds_indexed += s.seeds_indexed;
	}
};

/* Background thread of a worker in the parallel workflow that writes the current unix time
   as a 32 bit integer into the slot of the worker (offset 4 * worker_id) of the heartbeat
   file, once at startup and then every 15 minutes. Allows to detect stalled or killed workers:
   after each write, the slots of the other workers are compared against the current time and
   the workers more than 20 minutes behind are recorded in the Job (Job::stale_workers()).
   If a write fails, the thread terminates and check() rethrows the error in the calling thread.
   When the worker exits gracefully (the destructor does not run due to an exception), DONE is
   written into its slot so that it is not reported as stale by the other workers. */
struct Job;

struct Heartbeat {
	static const uint32_t DONE = 0xFFFFFFFF;
	Heartbeat(const std::string& file_name, Job& job);
	~Heartbeat();
	void check() const {
		if (failed_.load(std::memory_order_acquire))
			std::rethrow_exception(exception_);
	}
private:
	void run();
	void update_stale_workers(uint32_t now);
	FileStack file_;
	Job& job_;
	const int64_t worker_id_;
	std::mutex mtx_;
	std::condition_variable cv_;
	bool stop_;
	std::exception_ptr exception_;
	std::atomic<bool> failed_;
	const int uncaught_exceptions_;
	std::thread thread_;
};

struct Job {

	Job() :
		mem_limit(Util::String::interpret_number(config.memory_limit.get(DEFAULT_MEMORY_LIMIT))),
		max_oid_(std::numeric_limits<OId>::max()),
		base_dir_(absolute_path(config.tmpdir).first),
		round_(0),
		start_(std::chrono::system_clock::now()),
		steps_finished_(0)
	{
		register_temp_dir(base_dir_);
		make_temp_dir(base_dir());
		log_file_.reset(new FileStack(base_dir_ + PATH_SEPARATOR + "diamond_job.log", *this));
		// In recover mode, the job only inspects the temp directory, so it does not count as a worker.
		if (config.recover) {
			worker_id_ = -1;
			reset_workers();
		}
		else {
			Atomic worker_id(base_dir_ + PATH_SEPARATOR + "worker_id", *this);
			worker_id_ = worker_id.fetch_add();
		}
		if (!config.parallel_tmpdir.empty() && !config.recover) {
			const std::string heartbeat_file = base_dir_ + PATH_SEPARATOR + "heartbeat";
			register_sync_file(heartbeat_file);
			heartbeat_.reset(new Heartbeat(heartbeat_file, *this));
		}
	}

	void finish() {
		heartbeat_.reset();
		log_file_.reset();
		for (const auto& file : sync_files_) {
			remove_tmp_file(file);
		}
		for (auto it = temp_dirs_.rbegin(); it != temp_dirs_.rend(); ++it) {
			rmdir(*it);
		}
		//rmdir(base_dir_);
	}

	int64_t worker_id() const {
		return worker_id_;
	}

	std::string root_dir() const {
		return base_dir_ + PATH_SEPARATOR;
	}

	std::string base_dir(int round = -1) const {
		return base_dir_ + PATH_SEPARATOR + "round" + std::to_string(round == -1 ? round_ : round) + PATH_SEPARATOR;
	}

	void log(const char* format, ...);
	void log_raw(const std::string& message);
	std::string log_prefix() const;
	void log(const ClusterStats& stats);

	void next_round() {
		++round_;
		make_temp_dir(base_dir());
	}

	int round() const {
		return round_;
	}

	void set_round_count(int n, const std::vector<std::string>& steps) {
		round_count_ = n;
		steps_ = steps;
	}

	int round_count() const {
		return round_count_;
	}

	bool last_round() const {
		return round_ == round_count_ - 1;
	}

	const std::vector<std::string>& steps() const {
		return steps_;
	}
	
	const std::string& current_round() const {
		return steps_[round_];
	}

	bool is_linear_round() const {
		return ends_with(current_round(), "_lin");
	}

	ClusterStats& stats() {
		return stats_;
	}

	OId max_oid() const {
		if (max_oid_ == std::numeric_limits<OId>::max())
			throw std::runtime_error("max_oid not set");
		return max_oid_;
	}

	void set_max_oid(OId max_oid) {
		max_oid_ = max_oid;
	}

	/* Maximum length of the sequence accessions in the input, determined when the accession
	   files are written. Serves as an upper bound for sizing the output tables. */
	size_t max_acc_len() const {
		if (max_acc_len_ == std::numeric_limits<size_t>::max())
			throw std::runtime_error("max_acc_len not set");
		return max_acc_len_;
	}

	void set_max_acc_len(size_t len) {
		max_acc_len_ = len;
	}

	void register_sync_file(const std::string& file_name) {
		sync_files_.push_back(file_name);
	}

	void register_temp_dir(const std::string& dir_name) {
		temp_dirs_.push_back(dir_name);
	}

	void make_temp_dir(const std::string& dir_name) {
		mkdir(dir_name);
		register_temp_dir(dir_name);
	}

	void finish_step() {
		++steps_finished_;
	}

	void check_heartbeat() const {
		if (heartbeat_)
			heartbeat_->check();
	}

	/* Ids of the other workers whose heartbeat was more than 20 minutes behind the current
	   time at the last heartbeat of this worker. */
	std::vector<int64_t> stale_workers() const {
		std::lock_guard<std::mutex> lock(stale_workers_mtx_);
		return stale_workers_;
	}

	void set_stale_workers(std::vector<int64_t>&& workers) {
		std::lock_guard<std::mutex> lock(stale_workers_mtx_);
		stale_workers_ = std::move(workers);
	}

	std::string stop_file() const {
		return root_dir() + "stop";
	}

	/* Returns true if the parallel computation has to be stopped, either because this worker
	   detected stale workers (it then writes the stop file into the temp directory) or because
	   the stop file exists. Once true, it stays true, so that all callers along the way out
	   see the same decision even if a stale worker writes a heartbeat again in the meantime. */
	bool stopped();

	/* Recover mode: removes the stop file and the heartbeat file and resets the worker id
	   counter to 0, so that the workers of the next run start like those of a new run. */
	void reset_workers();

	/* Makes stopped() return true from now on, without writing the stop file, so that the
	   callers return as for a stopped job. Ends the run in recover mode. */
	void stop() {
		stopped_ = true;
	}

	bool goon() {
		if (stopped())
			return false;
		const bool r = !config.single_step || steps_finished_ == 0;
		if (!r)
			log("Job step finished, terminating...");
		return r;
	}

	const uint64_t mem_limit;

private:

	OId max_oid_;
	size_t max_acc_len_ = std::numeric_limits<size_t>::max();
	std::string base_dir_;
	int64_t worker_id_;
	int round_, round_count_;
	std::vector<std::string> steps_;
	std::unique_ptr<FileStack> log_file_;
	std::chrono::system_clock::time_point start_;
	ClusterStats stats_;
	std::vector<std::string> sync_files_;
	std::vector<std::string> temp_dirs_;
	int steps_finished_;
	mutable std::mutex stale_workers_mtx_;
	std::vector<int64_t> stale_workers_;
	bool stopped_ = false;
	// Declared last so that the thread is stopped before the other members are destroyed.
	std::unique_ptr<Heartbeat> heartbeat_;

};

/* Recover mode (diamond recover): the code of the parallel workflow is run through, checking the
   state of each synchronized step in the temp directory instead of executing it. A completed step
   is passed through like a worker that waited for it. At the first step that has not been
   completed, the job is stopped. If that step has been begun, its lock and queue are reset to 0,
   so that the next run executes it again (the queue skips the items already finished). */

/* Checks a step executed by the worker that takes the lock and marked as completed by done.
   queue_file, if given, is the queue the step uses for its work items. Returns true if the step
   has been completed, otherwise resets lock and queue and stops the job. */
bool recover_step(Job& job, const char* step, const std::string& lock_file, const std::string& done_file, const std::string& queue_file = std::string());
/* Checks a step whose n work items are distributed by a queue and recorded in the finished list
   when done. Returns true if all items are finished, otherwise resets the queue and stops the job. */
bool recover_queue(Job& job, const char* step, const std::string& queue_file, const std::string& finished_file, int64_t n);
// Returns true if the finished list holds all n work items. Does not create the file.
bool queue_finished(Job& job, const std::string& finished_file, int64_t n);

std::pair<std::string, uint64_t> get_reps(Job& job, const std::string& round_minichunks, std::unique_ptr<std::vector<OId>> clustering);
void merge(Job& job, const VolumedFile& volumes, const VolumedFile& input_volumes, Header hdr_format);
std::vector<int> make_blocks(Job& job, VolumedFile& volumes, std::vector<std::unique_ptr<std::ofstream>>& out, std::vector<std::unique_ptr<std::ofstream>>& acc_out,
	const std::string& seqs_vf, const std::string& accs_vf);
std::pair<std::string, std::string> len_sort(Job& job, VolumedFile& volumes);
std::vector<OId> build_merged(Job& job);
/* Loads the complete member superblock of a linear round. It is kept in memory while the
   superblock is searched against all rep superblocks and used by the searches without
   copying, instead of reading it from disk for each block combination. */
std::shared_ptr<Block> load_member_block(Job& job, const VolumedFile& volumes, int64_t member_superblock);
/* member_block is the loaded superblock r for linear rounds, and ignored otherwise. */
void run_search(Job& job, const VolumedFile& volumes, int64_t r, int64_t i, const std::shared_ptr<Block>& member_block, std::string base_dir, std::unique_ptr<std::vector<BitVector>>& seed_hit_table);
void configure_round(Job& job, uint64_t letter_count);
/* Masking algorithm applied to the input sequences when the length sorted minichunks are
   written. Everything downstream (seed counting, seed index building, alignment) reads the
   minichunks and thus works on already masked sequences. */
MaskingAlgo input_masking_algo();
bool use_lin_index(const Job& job);
std::string lin_index_file(const std::string& volume_path);
void build_lin_indices(Job& job, const VolumedFile& volumes);
/* Text file written next to the volume list of the superblocks, holding one line per
   superblock: the path of the superblock followed by the estimated number of distinct
   seeds of its most frequent seed shape, i.e. by the number of pivots that its seed index
   will hold. The hash tables of the index are sized from this count. */
std::string seed_count_file(const std::string& volume_list_file);
std::unordered_map<std::string, uint64_t> read_seed_counts(const std::string& file_name);
void remove_lin_indices(const VolumedFile& volumes);
std::pmr::unordered_map<OId, std::pmr::string> read_mapping_table(Job& job, const Volume& vol, size_t v, std::pmr::memory_resource& pool, bool remove);
std::string make_merged_blocks(Job& job, const std::string& minichunks, const std::string& merged_dir, uint64_t letter_count);
std::pair<std::vector<std::vector<HyperLogLog>>, std::vector<uint64_t>> count_distinct_seeds(Job& job, const VolumedFile& minichunks, const uint64_t letter_count);

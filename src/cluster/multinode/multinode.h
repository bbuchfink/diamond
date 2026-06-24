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
#include <chrono>
#include "basic/config.h"
#include "util/string/string.h"
#include "util/parallel/filestack.h"
#include "util/parallel/atomic.h"
#include "util/system/system.h"
#include "volume.h"

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

struct Job {

	Job() :
		mem_limit(Util::String::interpret_number(config.memory_limit.get(DEFAULT_MEMORY_LIMIT))),
		max_oid_(std::numeric_limits<OId>::max()),
		base_dir_(config.tmpdir),
		round_(0),
		start_(std::chrono::system_clock::now())
	{
		register_temp_dir(base_dir_);
		make_temp_dir(base_dir());
		log_file_.reset(new FileStack(base_dir_ + PATH_SEPARATOR + "diamond_job.log", *this));
		Atomic worker_id(base_dir_ + PATH_SEPARATOR + "worker_id", *this);
		worker_id_ = worker_id.fetch_add();
	}

	void finish() {
		//log("Cleaning up");
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
	void log(const ClusterStats& stats);

	void next_round() {
		++round_;
		make_temp_dir(base_dir());
	}

	int round() const {
		return round_;
	}

	void set_round(int64_t input_count) {
		input_count_.push_back(input_count);
	}

	uint64_t sparse_input_count(int round) const {
		return input_count_[round];
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

	void register_sync_file(const std::string& file_name) {
		/*if (!ends_with(file_name, "diamond_job.log"))
			log("Temp file: %s", file_name.c_str());
		else
			return;*/
		sync_files_.push_back(file_name);
	}

	void register_temp_dir(const std::string& dir_name) {
		temp_dirs_.push_back(dir_name);
	}

	void make_temp_dir(const std::string& dir_name) {
		mkdir(dir_name);
		register_temp_dir(dir_name);
	}

	const uint64_t mem_limit;

private:

	OId max_oid_;
	std::string base_dir_;
	int64_t worker_id_;
	int round_, round_count_;
	std::vector<std::string> steps_;
	std::unique_ptr<FileStack> log_file_;
	std::chrono::system_clock::time_point start_;
	std::vector<uint64_t> input_count_;
	ClusterStats stats_;
	std::vector<std::string> sync_files_;
	std::vector<std::string> temp_dirs_;

};

std::pair<std::string, uint64_t> get_reps(Job& job, const VolumedFile& volumes);
void write_representatives(Job& job, const VolumedFile& volumes, const std::vector<OId>& merged);
void merge(Job& job, const VolumedFile& volumes, Header hdr_format, const std::vector<OId>& merged);
//void extend(Job& job, std::vector<std::pair<OId, OId>>& out, const VolumedFile& volumes);
std::string len_sort(Job& job, VolumedFile& volumes);
std::vector<OId> build_merged(Job& job);
void run_search(Job& job, const VolumedFile& volumes, int64_t r, int64_t i, std::string base_dir, std::unique_ptr<std::vector<BitVector>>& seed_hit_table);
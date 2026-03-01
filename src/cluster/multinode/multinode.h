/****
Copyright (C) 2012-2026 Benjamin J. Buchfink <buchfink@gmail.com>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

	http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
****/
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <chrono>
#include "basic/config.h"
#include "util/string/string.h"
#include "util/parallel/filestack.h"
#include "basic/const.h"
#include "util/parallel/atomic.h"
#include "masking/masking.h"
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

	Job(OId max_oid, size_t volumes) :
		max_oid(max_oid),
		volumes(volumes),
		mem_limit(Util::String::interpret_number(config.memory_limit.get(DEFAULT_MEMORY_LIMIT))),
		base_dir_(config.parallel_tmpdir + PATH_SEPARATOR + "diamond-tmp-" + Const::version_string),
		round_(0),
		start_(std::chrono::system_clock::now())
	{
		mkdir(base_dir_);
		mkdir(base_dir());
		log_file_.reset(new FileStack(base_dir_ + PATH_SEPARATOR + "diamond_job.log"));
		Atomic worker_id(base_dir_ + PATH_SEPARATOR + "worker_id");
		worker_id_ = worker_id.fetch_add();
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
		mkdir(base_dir());
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

	void set_round_count(int n) {
		round_count_ = n;
	}

	int round_count() const {
		return round_count_;
	}

	bool last_round() const {
		return round_ == round_count_ - 1;
	}

	const OId max_oid;
	const size_t volumes;
	const uint64_t mem_limit;

private:

	std::string base_dir_;
	int64_t worker_id_;
	int round_, round_count_;
	std::unique_ptr<FileStack> log_file_;
	std::chrono::system_clock::time_point start_;
	std::vector<uint64_t> input_count_;

};

std::string get_reps(Job& job, const VolumedFile& volumes);
void merge(Job& job, const VolumedFile& volumes);
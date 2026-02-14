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
#include <vector>
#include <string>
#include <mutex>
#include "util/system/system.h"
#include "util/tsv/file.h"
#include "util/parallel/filestack.h"
#include "util/parallel/atomic.h"
#include "util/io/compressed_buffer.h"
#include "util/io/input_file.h"
#include "util/io/serializer.h"
#include "data/block/block.h"
#include "volume.h"
#include "util/algo/hash.h"
#include "file_array.h"

#ifdef WIN32
#else
#include <sys/stat.h>
#endif

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

#pragma pack(1)
struct PairEntry {
	static constexpr bool POD = true;
	PairEntry() :
		rep_oid(),
		member_oid(),
		rep_len(),
		member_len()
	{}
	PairEntry(int64_t rep_oid, int64_t member_oid, int32_t rep_len, int32_t member_len) :
		rep_oid(rep_oid),
		member_oid(member_oid),
		rep_len(rep_len),
		member_len(member_len)
	{}
	uint64_t key() const {
		return hash64(rep_oid);
	}
	bool operator<(const PairEntry& e) const {
		return rep_oid < e.rep_oid || (rep_oid == e.rep_oid && member_oid < e.member_oid);
	}
	friend void serialize(const PairEntry& e, CompressedBuffer& buf) {
		buf.write(e.rep_oid);
		buf.write(e.member_oid);
		buf.write(e.rep_len);
		buf.write(e.member_len);
	}
	friend void serialize(const PairEntry& e, Serializer& out) {
		out.write(e.rep_oid);
		out.write(e.member_oid);
		out.write(e.rep_len);
		out.write(e.member_len);
	}
	friend void deserialize(InputFile& in, PairEntry& e) {
		in.read(&e.rep_oid);
		in.read(&e.member_oid);
		in.read(&e.rep_len);
		in.read(&e.member_len);
	}
	struct Key {
		uint64_t operator()(const PairEntry& e) const {
			return e.rep_oid;
		}
	};
	uint64_t rep_oid, member_oid;
	uint32_t rep_len, member_len;
} PACKED_ATTRIBUTE;

struct PairEntryShort {
	static constexpr bool POD = true;
	PairEntryShort() :
		rep_oid(),
		member_oid()
	{}
	PairEntryShort(uint64_t rep_oid, uint64_t member_oid) :
		rep_oid(rep_oid),
		member_oid(member_oid)
	{}
	friend void serialize(const PairEntryShort& e, CompressedBuffer& buf) {
		buf.write(e.rep_oid);
		buf.write(e.member_oid);
	}
	friend void serialize(const PairEntryShort& e, Serializer& out) {
		out.write(e.rep_oid);
		out.write(e.member_oid);
	}
	friend void deserialize(InputFile& in, PairEntryShort& e) {
		in.read(&e.rep_oid);
		in.read(&e.member_oid);
	}
	struct Key {
		uint64_t operator()(const PairEntryShort& e) const {
			return e.rep_oid;
		}
	};
	uint64_t rep_oid, member_oid;
} PACKED_ATTRIBUTE;

struct Edge {
	static constexpr bool POD = true;
	Edge(uint64_t rep_oid, uint64_t member_oid, int32_t rep_len, int32_t member_len) :
		rep_oid(rep_oid),
		member_oid(member_oid),
		rep_len(rep_len),
		member_len(member_len)
	{}
	Edge() :
		rep_oid(),
		member_oid(),
		rep_len(),
		member_len()
	{}
	uint64_t key() const {
		return hash64(member_oid);
	}
	bool operator<(const Edge& e) const {
		return member_oid < e.member_oid || (member_oid == e.member_oid && (rep_len > e.rep_len || (rep_len == e.rep_len && rep_oid < e.rep_oid)));
	}
	friend void serialize(const Edge& e, CompressedBuffer& buf) {
		buf.write(e.rep_oid);
		buf.write(e.member_oid);
		buf.write(e.rep_len);
		buf.write(e.member_len);
	}
	friend void serialize(const Edge& e, Serializer& out) {
		out.write(e.rep_oid);
		out.write(e.member_oid);
		out.write(e.rep_len);
		out.write(e.member_len);
	}
	friend void deserialize(InputFile& in, Edge& e) {
		in.read(&e.rep_oid);
		in.read(&e.member_oid);
		in.read(&e.rep_len);
		in.read(&e.member_len);
	}
	struct Member {
		uint64_t operator()(const Edge& e) const {
			return e.member_oid;
		}
	};
	uint64_t rep_oid, member_oid;
	uint32_t rep_len, member_len;
} PACKED_ATTRIBUTE;

struct Assignment {
	static constexpr bool POD = true;
	Assignment() :
		member_oid(),
		rep_oid()
	{}
	Assignment(uint64_t member_oid, uint64_t rep_oid) :
		member_oid(member_oid),
		rep_oid(rep_oid)
	{}
	friend void serialize(const Assignment& e, CompressedBuffer& buf) {
		buf.write(e.member_oid);
		buf.write(e.rep_oid);
	}
	friend void serialize(const Assignment& e, Serializer& out) {
		out.write(e.member_oid);
		out.write(e.rep_oid);
	}
	friend void deserialize(InputFile& in, Assignment& e) {
		in.read(&e.member_oid);
		in.read(&e.rep_oid);
	}
	uint64_t member_oid, rep_oid;
} PACKED_ATTRIBUTE;
#pragma pack()

inline std::string base_path(const std::string& file_name) {
	static const std::string s = std::string(&PATH_SEPARATOR, 1) + "0" + PATH_SEPARATOR + "bucket.tsv";
	if (file_name.compare(file_name.length() - s.length(), s.length(), s) != 0)
		throw std::runtime_error("base_path");
	return file_name.substr(0, file_name.length() - s.length());
}

struct Job {

	Job(OId max_oid, size_t volumes):
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
		return base_dir_ + PATH_SEPARATOR + "round" + std::to_string(round == -1 ? round_ : round);
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

RadixedTable align(Job& job, int chunk_count, int64_t db_size);
std::string cluster(Job& job, const RadixedTable& edges, const VolumedFile& volumes);
std::string cluster_bidirectional(Job& job, const RadixedTable& edges, const VolumedFile& volumes);
void output(Job& job, const VolumedFile& volumes);
RadixedTable build_seed_table(Job& job, const VolumedFile& volumes, int shape);
RadixedTable build_pair_table(Job& job, const RadixedTable& seed_table, int shape, int64_t max_oid, FileArray& output_files);
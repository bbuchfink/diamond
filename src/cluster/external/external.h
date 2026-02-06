/****
Copyright � 2012-2026 Benjamin J. Buchfink <buchfink@gmail.com>

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
#include <atomic>
#include <stdarg.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <mutex>
#include "util/system/system.h"
#include "util/io/output_file.h"
#include "util/tsv/file.h"
#include "util/algo/partition.h"
#include "util/parallel/filestack.h"
#include "util/parallel/atomic.h"
#include "util/io/compressed_buffer.h"
#include "util/io/input_file.h"
#include "util/io/serializer.h"
#include "util/parallel/simple_thread_pool.h"
#include "data/block/block.h"

#ifdef WIN32
#else
#include <sys/stat.h>
#endif

const uint64_t RADIX_BITS = 8;
const int_fast16_t RADIX_COUNT = INT64_C(1) << RADIX_BITS;
const int_fast64_t MAX_FILE_SIZE = 1 * 1024 * 1024 * 1024;

#pragma pack(1)
struct PairEntry {
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
	int64_t key() const {
		return rep_oid;
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
		int64_t operator()(const PairEntry& e) const {
			return e.rep_oid;
		}
	};
	int64_t rep_oid, member_oid;
	int32_t rep_len, member_len;
} PACKED_ATTRIBUTE;

struct PairEntryShort {
	PairEntryShort() :
		rep_oid(),
		member_oid()
	{}
	PairEntryShort(int64_t rep_oid, int64_t member_oid) :
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
		int64_t operator()(const PairEntryShort& e) const {
			return e.rep_oid;
		}
	};
	int64_t rep_oid, member_oid;
} PACKED_ATTRIBUTE;

struct Edge {
	Edge(int64_t rep_oid, int64_t member_oid, int32_t rep_len, int32_t member_len) :
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
	int64_t key() const {
		return member_oid;
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
		int64_t operator()(const Edge& e) const {
			return e.member_oid;
		}
	};
	int64_t rep_oid, member_oid;
	int32_t rep_len, member_len;
} PACKED_ATTRIBUTE;

struct Assignment {
	Assignment() :
		member_oid(),
		rep_oid()
	{}
	Assignment(int64_t member_oid, int64_t rep_oid) :
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
	int64_t member_oid, rep_oid;
} PACKED_ATTRIBUTE;
#pragma pack()

inline std::string path(const std::string& file_name) {
	return file_name.substr(0, file_name.find_last_of(PATH_SEPARATOR));
}

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

	std::string base_dir(int round = -1) const {
		return base_dir_ + PATH_SEPARATOR + "round" + std::to_string(round == -1 ? round_ : round);
	}

	void log(const char* format, ...);

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

private:
	
	std::string base_dir_;
	int64_t worker_id_;
	int round_, round_count_;
	std::unique_ptr<FileStack> log_file_;
	std::chrono::system_clock::time_point start_;
	std::vector<uint64_t> input_count_;
	
};

struct Volume {
	Volume() :
		path(),
		oid_begin(0),
		oid_end(0),
		record_count(0)
	{}
	Volume(const std::string& path, OId oid_begin, OId oid_end, OId record_count) :
		path(path),
		oid_begin(oid_begin),
		oid_end(oid_end),
		record_count(record_count)
	{}
	std::string path;
	OId oid_begin, oid_end, record_count;
	bool operator<(size_t oid) const {
		return oid_end <= oid;
	}
	OId oid_range() const {
		return oid_end - oid_begin;
	}
	bool operator<(const Volume& v) const {
		return oid_begin < v.oid_begin;
	}
	friend std::istream& operator>>(std::istream& str, Volume& v) {
		std::string line;
		v.oid_begin = v.oid_end = 0;
		if (!std::getline(str, line)) return str;		
		std::istringstream row(line);
		row >> v.path;
		if(!row)
			throw std::runtime_error("Format error in VolumedFile");
		row >> v.record_count;
		if (!row)
			throw std::runtime_error("Format error in VolumedFile");
		row >> v.oid_begin >> v.oid_end;
		return str;
	}
};

struct VolumedFile : public std::vector<Volume> {
	VolumedFile(const std::string& file_name):
		list_file_(file_name),
		records_(0),
		max_oid_(0)
	{
		std::ifstream volume_file(file_name);
		int64_t oid = 0;
		Volume v;
		while(volume_file) {
			volume_file >> v;
			if (!volume_file)
				break;
			if (v.oid_begin == 0 && v.oid_end == 0) {
				v.oid_begin = oid;
				v.oid_end = oid + v.record_count;
			}
			push_back(v);
			oid += v.record_count;
			records_ += v.record_count;
			max_oid_ = std::max(max_oid_, v.oid_end - 1);
		}
		std::sort(begin(), end());
	}
	OId sparse_records() const {
		return records_;
	}
	OId max_oid() const {
		return max_oid_;
	}
	std::pair<std::vector<Volume>::const_iterator, std::vector<Volume>::const_iterator> find(OId oid_begin, OId oid_end) const {
		auto it = std::lower_bound(begin(), end(), oid_begin);
		if (it == end())
			throw std::runtime_error("OID out of bounds");
		auto end = it + 1;
		while (end < this->end() && end->oid_begin < oid_end)
			++end;
		return { it,end };
	}
	void remove() const {
		for (const Volume& v : *this)
			::remove(v.path.c_str());
		::remove(list_file_.c_str());
		rmdir(path(list_file_).c_str());
		
	}
private:
	const std::string list_file_;
	OId records_, max_oid_;
};

inline std::vector<std::string> read_list(const std::string& file_name) {
	std::vector<std::string> v;
	Util::Tsv::File file({ Util::Tsv::Type::STRING }, file_name);
	const Util::Tsv::Table table = file.read(config.threads_);
	v.reserve(table.size());
	for (int64_t i = 0; i < table.size(); ++i) {
		v.push_back(table[i].get<std::string>(0));
	}
	return v;
}

template<typename T>
struct InputBuffer {

	InputBuffer(const VolumedFile& f, int parts = config.threads_):
		size_(f.sparse_records()),
		data_(new T[size_]),
		part_(size_, parts)
	{
		std::atomic<int64_t> next(0);
		SimpleThreadPool pool;
		auto worker = [&](const std::atomic<bool>& stop) {
			int64_t v;
			while (!stop.load(std::memory_order_relaxed) && (v = next.fetch_add(1, std::memory_order_relaxed), v < (int64_t)f.size())) {
				InputFile in(f[v].path);
				in.read(&data_[f[v].oid_begin], f[v].record_count);
				in.close();
			}
			};
		for (int i = 0; i < std::min(config.threads_, (int)f.size()); ++i)
			pool.spawn(worker);
		pool.join_all();
	}

	size_t size() const {
		return size_;
	}

	size_t byte_size() const {
		return size_ * sizeof(T);
	}

	T* begin() {
		return data_.get();
	}

	T* end() {
		return data_.get() + size_;
	}

	const T* end() const {
		return data_.get() + size_;
	}

	const T* begin(int part) const {
		const T* begin = data_.get() + part_.begin(part), * end = this->end();
		if (part > 0)
			while (begin < end && begin[-1].key() == begin[0].key())
				++begin;
		return begin;
	}

	const T* end(int part) const {
		const T* ptr = data_.get() + part_.end(part), *end = this->end();
		while (ptr < end && ptr[-1].key() == ptr[0].key())
			++ptr;
		return ptr;
	}

	int64_t parts() const {
		return part_.parts;
	}

	const T& front() const {
		return data_[0];
	}

	const T& back() const {
		return data_[size_ - 1];
	}

private:

	const size_t size_;
	std::unique_ptr<T[]> data_;
	const Partition<int64_t> part_;

};

std::vector<std::string> align(Job& job, int chunk_count, int64_t db_size);
std::string cluster(Job& job, const std::vector<std::string>& edges, const VolumedFile& volumes);
std::string cluster_bidirectional(Job& job, const std::vector<std::string>& edges, const VolumedFile& volumes);
void output(Job& job);
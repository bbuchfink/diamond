/****
Copyright © 2012-2026 Benjamin J. Buchfink <buchfink@gmail.com>

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
#include "data/block/block.h"

#ifdef WIN32
#else
#include <sys/stat.h>
#endif

const uint64_t RADIX_BITS = 8;
const int_fast16_t RADIX_COUNT = INT64_C(1) << RADIX_BITS;
const int_fast64_t MAX_FILE_SIZE = 1 * 1024 * 1024 * 1024;

struct PairEntry {
	PairEntry() :
		rep_oid(),
		member_oid(),
		rep_len(),
		member_len()
	{}
	PairEntry(int_fast64_t rep_oid, int_fast64_t member_oid, int_fast32_t rep_len, int_fast32_t member_len) :
		rep_oid(rep_oid),
		member_oid(member_oid),
		rep_len(rep_len),
		member_len(member_len)
	{}
	int_fast64_t key() const {
		return rep_oid;
	}
	bool operator<(const PairEntry& e) const {
		return rep_oid < e.rep_oid || (rep_oid == e.rep_oid && member_oid < e.member_oid);
	}
	struct Key {
		int64_t operator()(const PairEntry& e) const {
			return e.rep_oid;
		}
	};
	int_fast64_t rep_oid, member_oid;
	int_fast32_t rep_len, member_len;
};

struct PairEntryShort {
	PairEntryShort() :
		rep_oid(),
		member_oid()
	{}
	PairEntryShort(int64_t rep_oid, int64_t member_oid) :
		rep_oid(rep_oid),
		member_oid(member_oid)
	{}
	struct Key {
		int64_t operator()(const PairEntryShort& e) const {
			return e.rep_oid;
		}
	};
	int64_t rep_oid, member_oid;
};

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
	struct Member {
		int64_t operator()(const Edge& e) const {
			return e.member_oid;
		}
	};
	int64_t rep_oid, member_oid;
	int32_t rep_len, member_len;
};

struct Assignment {
	int64_t member_oid, rep_oid;
};

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

struct FileArray {

	FileArray(const std::string& base_dir, int size, int64_t worker_id, int64_t max_file_size = MAX_FILE_SIZE) :
		max_file_size(max_file_size),
		size_(size),
		worker_id_(worker_id),
		base_dir(base_dir),
		mtx_(size),
		records_(size, 0),
		bytes_(size, 0),
		next_(size, 1)
	{
		for (int64_t i = 0; i < size; ++i) {
			const std::string dir = base_dir + PATH_SEPARATOR + std::to_string(i) + PATH_SEPARATOR;
			mkdir(dir);
			output_files_.push_back(new OutputFile(dir + "worker_" + std::to_string(worker_id) + "_volume_0"));
			bucket_files_.emplace_back(new FileStack(dir + "bucket.tsv"));
		}
	}

	~FileArray() {
		for (int64_t i = 0; i < size_; ++i) {
			output_files_[i]->close();
			if (records_[i] > 0)
				bucket_files_[i]->push(output_files_[i]->file_name() + '\t' + std::to_string(records_[i]));
			else
				::remove(output_files_[i]->file_name().c_str());
			delete output_files_[i];
		}
	}

	bool write(int i, const char* ptr, int64_t count, int64_t records) {
		std::lock_guard<std::mutex> lock(mtx_[i]);
		output_files_[i]->write(ptr, count);
		records_[i] += records;
		bytes_[i] += count;
		if (bytes_[i] >= max_file_size) {
			bucket_files_[i]->push(output_files_[i]->file_name() + '\t' + std::to_string(records_[i]));
			records_[i] = 0;
			bytes_[i] = 0;
			output_files_[i]->close();
			delete output_files_[i];
			output_files_[i] = new OutputFile(base_dir + PATH_SEPARATOR + std::to_string(i) + PATH_SEPARATOR + "worker_" + std::to_string(worker_id_) + "_volume_" + std::to_string(next_[i]++));
			return true;
		}
		return false;
	}

	int64_t records(int i) const {
		return records_[i];
	}

	std::string bucket(int i) const {
		return bucket_files_[i]->file_name();
	}
	
	std::vector<std::string> buckets() const {
		std::vector<std::string> buckets;
		buckets.reserve(size_);
		for (int i = 0; i < size_; ++i)
			buckets.push_back(bucket(i));
		return buckets;
	}

	std::string file_name(int i) {
		return output_files_[i]->file_name();
	}

private:

	const int64_t max_file_size;
	const int size_;
	const int64_t worker_id_;
	const std::string base_dir;
	std::vector<OutputFile*> output_files_;
	std::vector<std::mutex> mtx_;
	std::vector<int64_t> records_, bytes_, next_;
	std::vector<std::unique_ptr<FileStack>> bucket_files_;

};

/*template<typename T, int N>
struct BufferArray {
	static constexpr int64_t BUF_SIZE = 4096;
	BufferArray(FileArray& file_array) :
		file_array_(file_array)
	{}
	void write(int radix, const T& x) {
		data_[radix].push_back(x);
		if (data_[radix].size() >= BUF_SIZE) {
			file_array_.write(radix, data_[radix].data(), data_[radix].size());
			data_[radix].clear();
		}
	}
	~BufferArray() {
		for (int i = 0; i < N; ++i)
			file_array_.write(i, data_[i].data(), data_[i].size());
	}
private:
	std::array<std::vector<T>, N> data_;
	FileArray& file_array_;
};*/


struct BufferArray {
	static constexpr int64_t BUF_SIZE = 65536;
	BufferArray(FileArray& file_array, int size) :
		data_(size),
		records_(size, 0),
		file_array_(file_array)
	{
	}
	template<typename T>
	bool write(int radix, const T* ptr, int64_t n, int64_t record_count = -1) {
		data_[radix].write((const char*)ptr, n * sizeof(T));
		records_[radix] += record_count >= 0 ? record_count : n;
		if (data_[radix].size() >= BUF_SIZE) {
			data_[radix].finish();
			const bool r = file_array_.write(radix, data_[radix].data(), data_[radix].size(), records_[radix]);
			data_[radix].clear();
			records_[radix] = 0;
			return r;
		}
		return false;
	}
	template<typename T>
	bool write(int radix, const T& x) {
		return write(radix, &x, 1);
	}
	~BufferArray() {
		for (int i = 0; i < (int)data_.size(); ++i) {
			data_[i].finish();
			file_array_.write(i, data_[i].data(), data_[i].size(), records_[i]);
		}
	}
private:
	std::vector<CompressedBuffer> data_;
	std::vector<int64_t> records_;
	FileArray& file_array_;
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
		auto worker = [&]() {
			int64_t v;
			while (v = next.fetch_add(1, std::memory_order_relaxed), v < (int64_t)f.size()) {
				InputFile in(f[v].path);
				in.read(&data_[f[v].oid_begin], f[v].record_count);
				in.close();
			}
			};
		std::vector<std::thread> t;
		for (int i = 0; i < std::min(config.threads_, (int)f.size()); ++i)
			t.emplace_back(worker);
		for (auto& i : t)
			i.join();
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
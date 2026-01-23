/****
Copyright © 2013-2025 Benjamin J. Buchfink

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

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
#include "basic/sequence.h"
#include "data/block/block.h"
#include "data/sequence_file.h"

#ifdef WIN32
#else
#include <sys/stat.h>
#endif

namespace Cluster {

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

	Job():
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

	void set_round(size_t volumes, int64_t input_count) {
		volume_count_.push_back(volumes);
		input_count_.push_back(input_count);
	}

	size_t volume_count(int round) const {
		return volume_count_[round];
	}

	uint64_t input_count(int round) const {
		return input_count_[round];
	}

	void set_round_count(int n) {
		round_count_ = n;
	}

	bool last_round() const {
		return round_ == round_count_ - 1;
	}

private:

	std::string base_dir_;
	int64_t worker_id_;
	int round_, round_count_;
	std::unique_ptr<FileStack> log_file_;
	std::chrono::system_clock::time_point start_;
	std::vector<size_t> volume_count_;
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
	Volume(const std::string& path, size_t oid_begin, size_t record_count) :
		path(path),
		oid_begin(oid_begin),
		record_count(record_count)
	{}
	std::string path;
	size_t oid_begin, record_count;
	size_t oid_end() const {
		return oid_begin + record_count;
	}
	bool operator<(size_t oid) const {
		return oid_end() <= oid;
	}
};

struct VolumedFile : public std::vector<Volume> {
	VolumedFile(const std::string& file_name):
		list_file_(file_name)
	{
		Util::Tsv::File volume_file({ Util::Tsv::Type::STRING, Util::Tsv::Type::INT64 }, file_name);
		const Util::Tsv::Table volume_table = volume_file.read(config.threads_);
		reserve(volume_table.size());
		int64_t oid = 0;
		for (int64_t i = 0; i < volume_table.size(); ++i) {
			const int64_t n = volume_table[i].get<int64_t>(1);
			emplace_back(volume_table[i].get<std::string>(0), oid, n);
			oid += n;
		}
	}
	size_t records() const {
		return empty() ? 0 : back().oid_end();
	}
	std::pair<std::vector<Volume>::const_iterator, std::vector<Volume>::const_iterator> find(int64_t oid_begin, OId oid_end) const {
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
		size_(f.records()),
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

struct ChunkSeqs {

	ChunkSeqs(const std::string& chunk_path):
		seq_file_(chunk_path + "bucket.tsv"),
		oid_count_(0),
		letter_count_(0)
	{
		std::atomic<int> next(0);
		seq_blocks_.resize(seq_file_.size());
		oid2seq_.resize(seq_file_.size());
		oid_range_.resize(seq_file_.size());
		//std::mutex mtx;
		auto worker = [&]() {
			int i;
			while (i = next.fetch_add(1, std::memory_order_relaxed), i < (int)seq_file_.size()) {
				std::unique_ptr<SequenceFile> in(SequenceFile::auto_create({ seq_file_[i].path }));
				in->flags() |= SequenceFile::Flags::SEQS | SequenceFile::Flags::TITLES;
				seq_blocks_[i] = in->load_seqs(INT64_MAX);
				in->close();
				Block& b = *seq_blocks_[i];
				const StringSet& ids = b.ids();
				const SequenceSet& seqs = b.seqs();
				const BlockId n = ids.size();
				std::unordered_map<int64_t, Sequence>& m = *(oid2seq_[i] = new std::unordered_map<int64_t, Sequence>());
				m.reserve(n);
				int64_t oid_min = INT64_MAX, oid_max = INT64_MIN;
				for (BlockId j = 0; j < n; ++j) {
					const int64_t oid = std::atoll(ids[j]);
					oid_min = std::min(oid_min, oid);
					oid_max = std::max(oid_max, oid);
					m[oid] = seqs[j];
				}
				{
					//std::lock_guard<std::mutex> lock(mtx);
					//range2block_.insert(std::make_pair(std::pair<int64_t, int64_t>(oid_min, oid_max), i));
					oid_range_[i] = { oid_min, oid_max };
				}
			}
			};
		std::vector<std::thread> workers;
		for (int i = 0; i < std::min(config.threads_, int(seq_file_.size())); ++i)
			workers.emplace_back(worker);
		for (auto& t : workers)
			t.join();
		int64_t oid_count = 0, letter_count = 0;
		for (const Block* b : seq_blocks_) {
			oid_count_ += b->seqs().size();
			letter_count_ += b->seqs().letters();
		}
	}

	~ChunkSeqs() {
		std::vector<std::thread> workers;
		std::atomic<int> next(0);
		auto worker = [&]() {
			int i;
			while(i = next.fetch_add(1, std::memory_order_relaxed), i < (int)seq_blocks_.size()) {
				delete seq_blocks_[i];
				delete oid2seq_[i];
			}
			};
		for (int i = 0; i < std::min(config.threads_, int(seq_blocks_.size())); ++i)
			workers.emplace_back(worker);
		for (auto& t : workers)
			t.join();
		seq_file_.remove();
	}

	size_t oids() const {
		return oid_count_;
	}

	size_t letters() const {
		return letter_count_;
	}

	size_t volumes() const {
		return seq_blocks_.size();
	}

	Sequence operator[](int64_t oid) const {
		//const auto r = range2block_.equal_range(std::make_pair(oid, oid));
		//for (auto i = r.first; i != r.second; ++i) {
		for (int i = 0; i < (int)oid_range_.size(); ++i) {
			if (oid < oid_range_[i].first || oid > oid_range_[i].second)
				continue;
			std::unordered_map<int64_t, Sequence>::const_iterator it;
			if ((it = oid2seq_[i]->find(oid)) != oid2seq_[i]->end())
				return it->second;
		}
		throw std::out_of_range("ChunkSeqs");
	}

private:

	/*struct CmpRange {
		bool operator()(const std::pair<int64_t, int64_t>& a, const std::pair<int64_t, int64_t>& b) const {
			return a.second < b.first;
		}
	};*/

	VolumedFile seq_file_;
	int64_t oid_count_, letter_count_;
	std::vector<Block*> seq_blocks_;
	//std::multimap<std::pair<int64_t, int64_t>, int, ChunkSeqs::CmpRange> range2block_;
	std::vector<std::pair<int64_t, int64_t>> oid_range_;
	std::vector<std::unordered_map<int64_t, Sequence>*> oid2seq_;

};

std::vector<std::string> align(Job& job, int chunk_count, int64_t db_size);
std::string cluster(Job& job, const std::vector<std::string>& edges, const VolumedFile& volumes);
std::string cluster_bidirectional(Job& job, const std::vector<std::string>& edges, const VolumedFile& volumes);
void output(Job& job);

}
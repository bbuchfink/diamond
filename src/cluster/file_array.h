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

#pragma once
#include <stdexcept>
#include <string>
#include "radixed_table.h"
#include "util/io/compressed_buffer.h"
#include "util/algo/degree_partition.h"

const int RADIX_BITS = 8;
const uint64_t RADIX_COUNT = UINT64_C(1) << RADIX_BITS;
const uint64_t MAX_FILE_SIZE = 1 * 1024 * 1024 * 1024;

struct FileArray {

	FileArray(const std::string& base_dir, int size, int64_t worker_id, bool exclusive, int64_t max_file_size = MAX_FILE_SIZE) :
		exclusive_(exclusive),
		max_file_size(max_file_size),
		size_(size),
		worker_id_(worker_id),
		base_dir(base_dir),
		mtx_(size),
		records_(size, 0),
		records_total_(size, 0),
		bytes_(size, 0),
		next_(size, 1)
	{
		for (int64_t i = 0; i < size; ++i) {
			const std::string dir = base_dir + PATH_SEPARATOR + std::to_string(i) + PATH_SEPARATOR;
			mkdir(dir);
			output_files_.push_back(new File(dir + "worker_" + std::to_string(worker_id) + "_volume_0", "wb"));
			bucket_files_.emplace_back(new FileStack(dir + "bucket.tsv"));
		}
	}

	void close() {
		if (output_files_.empty())
			return;
		for (uint64_t i = 0; i < size_; ++i) {
			output_files_[i]->close();
			if (records_[i] > 0)
				bucket_files_[i]->push(output_files_[i]->name() + '\t' + std::to_string(records_[i]));
			else
				::remove(output_files_[i]->name().c_str());
			bucket_files_[i]->close();
			delete output_files_[i];
		}
		output_files_.clear();
	}

	~FileArray() {
		close();
	}

	bool write(uint64_t i, const char* ptr, size_t count, int64_t records) {
		std::lock_guard<std::mutex> lock(mtx_[i]);
		output_files_[i]->write(ptr, count);
		records_[i] += records;
		if (exclusive_)
			records_total_[i] += records;
		bytes_[i] += count;
		if (bytes_[i] >= max_file_size) {
			bucket_files_[i]->push(output_files_[i]->name() + '\t' + std::to_string(records_[i]));
			records_[i] = 0;
			bytes_[i] = 0;
			output_files_[i]->close();
			delete output_files_[i];
			output_files_[i] = new File(base_dir + PATH_SEPARATOR + std::to_string(i) + PATH_SEPARATOR + "worker_" + std::to_string(worker_id_) + "_volume_" + std::to_string(next_[i]++), "wb");
			return true;
		}
		return false;
	}

	uint64_t records(uint64_t i) const {
		return records_[i];
	}

	std::string bucket(uint64_t i) const {
		return bucket_files_[i]->file_name();
	}

	RadixedTable buckets(int shift) const {
		RadixedTable buckets(shift);
		buckets.reserve(size_);
		for (uint64_t i = 0; i < size_; ++i)
			buckets.emplace_back(bucket(i), exclusive_ ? records_total_[i] : Bucket::NIL, i << shift, (i + 1) << shift);
		return buckets;
	}

	RadixedTable buckets(const DegreePartition& p) {
		RadixedTable buckets(0);
		const auto b = p.buckets();
		if (b.size() != size_)
			throw std::runtime_error("FileArray::buckets");
		buckets.reserve(b.size());
		for (uint64_t i = 0; i < b.size(); ++i)
			buckets.emplace_back(bucket(i), exclusive_ ? records_total_[i] : Bucket::NIL, b[i].first_degree, b[i].last_degree + 1);
		return buckets;
	}

	std::string file_name(int i) {
		return output_files_[i]->name();
	}

	uint64_t records_total() const {
		if (!exclusive_)
			throw std::runtime_error("Total record count is only available in exclusive mode");
		uint64_t sum = 0;
		for (uint64_t i = 0; i < size_; ++i)
			sum += records_[i];
		return sum;
	}

private:

	const bool exclusive_;
	const int64_t max_file_size;
	const uint64_t size_;
	const int64_t worker_id_;
	const std::string base_dir;
	std::vector<File*> output_files_;
	std::vector<std::mutex> mtx_;
	std::vector<int64_t> records_, records_total_, bytes_, next_;
	std::vector<std::unique_ptr<FileStack>> bucket_files_;

};

struct BufferArray {

	static constexpr int64_t BUF_SIZE = 65536;

	BufferArray(FileArray& file_array, int size) :
		data_(size),
		records_(size, 0),
		file_array_(file_array)
	{
	}

	template<typename T>
	void write(uint64_t radix, const T* ptr, size_t n, int64_t record_count) {
		check_radix(radix);
		for (size_t i = 0; i < n; ++i)
			serialize(ptr[i], data_[radix]);
		records_[radix] += record_count;
		flush(radix);
	}

	void write(uint64_t radix, const char* ptr, size_t n) {
		check_radix(radix);
		data_[radix].write(ptr, n);
		records_[radix] += n;
		flush(radix);
	}

	void flush(uint64_t radix) {
		check_radix(radix);
		if (data_[radix].size() >= BUF_SIZE) {
			data_[radix].finish();
			file_array_.write(radix, data_[radix].data(), data_[radix].size(), records_[radix]);
			data_[radix].clear();
			records_[radix] = 0;
		}
	}

	template<typename T>
	void write(uint64_t radix, const T& x) {
		write(radix, &x, 1, 1);
	}

	template<typename T>
	void write_msb(const T& x) {
		const uint64_t key = x.key();
		const uint64_t radix = key >> (64 - RADIX_BITS);
		write(radix, &x, 1, 1);
	}

	void finish() {
		for (uint64_t i = 0; i < data_.size(); ++i) {
			data_[i].finish();
			file_array_.write(i, data_[i].data(), data_[i].size(), records_[i]);
		}
		data_.clear();
	}

	~BufferArray() {
		finish();
	}

private:

	void check_radix(uint64_t radix) const {
		if (radix >= (uint64_t)data_.size())
			throw std::out_of_range("BufferArray: radix out of range: " + std::to_string(radix)
				+ " (buckets: " + std::to_string(data_.size()) + ")");
	}

	std::vector<CompressedBuffer> data_;
	std::vector<int64_t> records_;
	FileArray& file_array_;

};
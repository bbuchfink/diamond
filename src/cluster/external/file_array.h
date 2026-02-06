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

	bool write(int i, const char* ptr, size_t count, int64_t records) {
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
	//bool write(int radix, const T* ptr, size_t n, int64_t record_count) {
	void write(int radix, const T * ptr, size_t n, int64_t record_count) {
		for (size_t i = 0; i < n; ++i)
			serialize(ptr[i], data_[radix]);
		records_[radix] += record_count;
		/*if (data_[radix].size() >= BUF_SIZE) {
			data_[radix].finish();
			const bool r = file_array_.write(radix, data_[radix].data(), data_[radix].size(), records_[radix]);
			data_[radix].clear();
			records_[radix] = 0;
			return r;
		}
		return false;*/
		flush(radix);
	}

	void write(int radix, const char* ptr, size_t n) {
		data_[radix].write(ptr, n);
		records_[radix] += n;
		flush(radix);
	}

	void flush(int radix) {
		if (data_[radix].size() >= BUF_SIZE) {
			data_[radix].finish();
			file_array_.write(radix, data_[radix].data(), data_[radix].size(), records_[radix]);
			data_[radix].clear();
			records_[radix] = 0;
		}
	}

	template<typename T>
	void write(int radix, const T& x) {
		write(radix, &x, 1, 1);
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
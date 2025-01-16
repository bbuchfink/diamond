/****
DIAMOND protein aligner
Copyright (C) 2013-2024 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen

Code developed by Benjamin Buchfink <buchfink@gmail.com>

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

#pragma once
#include <vector>
#include <assert.h>
#include <thread>
#include <iterator>
#include <atomic>
#include "io/input_file.h"
#include "log_stream.h"
#include "util/ptr_vector.h"
#include "io/async_file.h"
#include "io/input_stream_buffer.h"
#include "text_buffer.h"
#include "io/serialize.h"
#include "data_structures/writer.h"

template<typename T>
struct AsyncBuffer
{

	using Vector = std::vector<T>;
	using Key = typename SerializerTraits<T>::Key;
	static const int64_t ENTRY_SIZE = (int64_t)sizeof(T);

	AsyncBuffer(const std::vector<Key>& key_partition, const std::string &tmpdir, const SerializerTraits<T>& traits) :
		key_partition_(key_partition),
		traits_(traits),
		bins_processed_(0),
		total_disk_size_(0)
	{
		log_stream << "Async_buffer() " << key_partition.back() << std::endl;
		count_ = new std::atomic_size_t[key_partition.size()];
		for (size_t i = 0; i < key_partition.size(); ++i) {
			tmp_file_.push_back(new AsyncFile());
			count_[i] = (size_t)0;
		}
	}

	~AsyncBuffer() {
		delete[] count_;
	}

	Key begin(int bin) const
	{
		return bin == 0 ? 0 : key_partition_[bin - 1];
	}

	Key end(int bin) const
	{
		return key_partition_[bin];
	}

	int bins() const {
		return (int)key_partition_.size();
	}

	int bin(Key key) const {
		const int n = (int)key_partition_.size();
		for (int i = 0; i < n; ++i)
			if (key < key_partition_[i])
				return i;
		throw std::runtime_error("key_partition error");
	}

	struct Iterator : public Writer<T>
	{
		Iterator(AsyncBuffer& parent, size_t thread_num) :
			last_bin_(0),
			buffer_(parent.bins()),
			count_(parent.bins(), 0),
			parent_(parent)
		{
			ser_.reserve(parent.bins());
			for (int i = 0; i < parent.bins(); ++i) {
				ser_.emplace_back(buffer_[i], parent.traits_);
				out_.push_back(&parent.tmp_file_[i]);
			}
		}
		virtual Iterator& operator=(const T& x) override
		{
			if (SerializerTraits<T>::is_sentry(x)) {
				last_bin_ = parent_.bin(ser_.front().traits.key(x));
				if (buffer_[last_bin_].size() >= buffer_size)
					flush(last_bin_);
			}
			else
				++count_[last_bin_];
			assert(last_bin_ < parent_.bins());
			ser_[last_bin_] << x;
			return *this;
		}
		void flush(int bin)
		{
			out_[bin]->write(buffer_[bin].data(), buffer_[bin].size());
			buffer_[bin].clear();
		}
		virtual ~Iterator()
		{
			for (int bin = 0; bin < parent_.bins(); ++bin) {
				flush(bin);
				parent_.count_[bin] += count_[bin];
			}
		}
	private:
		enum { buffer_size = 65536 };
		int last_bin_;
		std::vector<TextBuffer> buffer_;
		std::vector<TypeSerializer<T>> ser_;
		std::vector<size_t> count_;
		std::vector<AsyncFile*> out_;
		AsyncBuffer &parent_;
	};

	void load(int64_t max_size) {
		max_size = std::max(max_size, (int64_t)1);
		auto worker = [this](int end) {
			for (; bins_processed_ < end; ++bins_processed_)
				load_bin(*data_next_, bins_processed_);
		};
		if (bins_processed_ == bins()) {
			data_next_ = nullptr;
			return;
		}
		int64_t size = count_[bins_processed_], current_size, disk_size = tmp_file_[bins_processed_].tell();
		int end = bins_processed_ + 1;
		while (end < bins() && (size + (current_size = count_[end])) * ENTRY_SIZE < max_size) {
			size += current_size;
			disk_size += tmp_file_[end].tell();
			++end;
		}
		log_stream << "Async_buffer.load() " << size << "(" << (double)size * sizeof(T) / (1 << 30) << " GB, " << (double)disk_size / (1 << 30) << " GB on disk)" << std::endl;
		total_disk_size_ += disk_size;
		data_next_ = new std::vector<T>;
		data_next_->reserve(size);
		input_range_next_.first = begin(bins_processed_);
		input_range_next_.second = this->end(end - 1);
		load_worker_ = new std::thread(worker, end);
	}

	std::tuple<std::vector<T>*, Key, Key> retrieve() {
		if (data_next_ != nullptr) {
			load_worker_->join();
			delete load_worker_;
		}
		return std::tuple<std::vector<T>*, Key, Key> { data_next_, input_range_next_.first, input_range_next_.second };
	}

	size_t total_disk_size() {
		return total_disk_size_;
	}

	int64_t bin_size(int i) const {
		return count_[i];
	}

private:

	void load_bin(std::vector<T> &out, size_t bin)
	{
		InputFile f(tmp_file_[bin], InputStreamBuffer::ASYNC);
		const size_t n = out.size();
		if (count_[bin] > 0) {
			auto it = std::back_inserter(out);
			TypeDeserializer<T>(f, traits_) >> it;
			if ((out.size() - n) != count_[bin])
				throw std::runtime_error("Mismatching hit count / possibly corrupted temporary file: " + f.file_name);
		}
		f.close_and_delete();
	}

	const std::vector<Key> key_partition_;
	const SerializerTraits<T> traits_;
	int bins_processed_;
	int64_t total_disk_size_;
	PtrVector<AsyncFile> tmp_file_;
	std::atomic_size_t *count_;
	std::pair<Key, Key> input_range_next_;
	std::vector<T>* data_next_;
	std::thread* load_worker_;

};
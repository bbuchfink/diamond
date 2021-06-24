/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen

Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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
#include <exception>
#include <assert.h>
#include <thread>
#include <tuple>
#include <iterator>
#include <atomic>
#include "io/temp_file.h"
#include "io/input_file.h"
#include "log_stream.h"
#include "../util/ptr_vector.h"
#include "io/async_file.h"
#include "io/input_stream_buffer.h"
#include "text_buffer.h"
#include "io/serialize.h"
#include "data_structures/writer.h"

template<typename T>
struct AsyncBuffer
{

	typedef std::vector<T> Vector;

	AsyncBuffer(size_t input_count, const std::string &tmpdir, unsigned bins, const SerializerTraits<T>& traits) :
		bins_(bins),
		bin_size_((input_count + bins_ - 1) / bins_),
		input_count_(input_count),
		traits_(traits),
		bins_processed_(0),
		total_disk_size_(0)
	{
		log_stream << "Async_buffer() " << input_count << ',' << bin_size_ << std::endl;
		count_ = new std::atomic_size_t[bins];
		for (unsigned i = 0; i < bins; ++i) {
			tmp_file_.push_back(new AsyncFile());
			count_[i] = (size_t)0;
		}
	}

	~AsyncBuffer() {
		delete[] count_;
	}

	size_t begin(size_t bin) const
	{
		return bin*bin_size_;
	}

	size_t end(size_t bin) const
	{
		return std::min((bin + 1)*bin_size_, input_count_);
	}

	struct Iterator : public Writer<T>
	{
		Iterator(AsyncBuffer &parent, size_t thread_num) :
			buffer_(parent.bins()),
			count_(parent.bins(), 0),
			parent_(parent)
		{
			ser_.reserve(parent.bins_);
			for (unsigned i = 0; i < parent.bins_; ++i) {
				ser_.emplace_back(buffer_[i], parent.traits_);
				out_.push_back(&parent.tmp_file_[i]);
			}
		}
		virtual Iterator& operator=(const T& x) override
		{
			const unsigned bin = ser_.front().traits.key(x) / parent_.bin_size_;
			if (SerializerTraits<T>::is_sentry(x)) {
				if (buffer_[bin].size() >= buffer_size)
					flush(bin);
			}
			else
				++count_[bin];
			assert(bin < parent_.bins());
			ser_[bin] << x;
			return *this;
		}
		void flush(unsigned bin)
		{
			out_[bin]->write(buffer_[bin].data(), buffer_[bin].size());
			buffer_[bin].clear();
		}
		virtual ~Iterator()
		{
			for (unsigned bin = 0; bin < parent_.bins_; ++bin) {
				flush(bin);
				parent_.count_[bin] += count_[bin];
			}
		}
	private:
		enum { buffer_size = 65536 };
		std::vector<TextBuffer> buffer_;
		std::vector<TypeSerializer<T>> ser_;
		std::vector<size_t> count_;
		std::vector<AsyncFile*> out_;
		AsyncBuffer &parent_;
	};

	void load(size_t max_size) {
		auto worker = [this](size_t end) {
			for (; bins_processed_ < end; ++bins_processed_)
				load_bin(*data_next_, bins_processed_);
		};
		if (bins_processed_ == bins_) {
			data_next_ = nullptr;
			return;
		}
		size_t size = count_[bins_processed_], end = bins_processed_ + 1, current_size, disk_size = tmp_file_[bins_processed_].tell();
		while (end < bins_ && (size + (current_size = count_[end])) * sizeof(T) < max_size) {
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

	std::tuple<std::vector<T>*, size_t, size_t> retrieve() {
		if (data_next_ != nullptr) {
			load_worker_->join();
			delete load_worker_;
		}
		return std::tuple<std::vector<T>*, size_t, size_t> { data_next_, input_range_next_.first, input_range_next_.second };
	}

	unsigned bins() const
	{
		return bins_;
	}

	size_t total_disk_size() {
		return total_disk_size_;
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

	const unsigned bins_;
	const size_t bin_size_, input_count_;
	const SerializerTraits<T> traits_;
	size_t bins_processed_, total_disk_size_;
	PtrVector<AsyncFile> tmp_file_;
	std::atomic_size_t *count_;
	std::pair<size_t, size_t> input_range_next_;
	std::vector<T>* data_next_;
	std::thread* load_worker_;

};
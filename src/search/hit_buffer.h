/****
DIAMOND protein aligner
Copyright (C) 2013-2025 Max Planck Society for the Advancement of Science e.V.
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
#include <assert.h>
#include <vector>
#include <thread>
#include <atomic>
#include <mutex>
#include "basic/config.h"
#include "util/io/mmap.h"
#include "util/text_buffer.h"
#include "hit.h"

namespace Search {

struct HitBuffer
{

	using Vector = std::vector<Hit>;
	using Key = uint32_t;
	static const int64_t ENTRY_SIZE = (int64_t)sizeof(Hit);

	HitBuffer(const std::vector<Key>& key_partition, const std::string& tmpdir, bool long_subject_offsets, int query_contexts);
	~HitBuffer() {
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

	struct Writer
	{
		Writer(HitBuffer& parent, size_t thread_num) :
			buffer_size(config.trace_pt_membuf ? 8192 : 65536),
			last_bin_(0),
			buffer_(config.trace_pt_membuf ? parent.bins() : 0),
			text_buffer_(config.trace_pt_membuf ? 0 : parent.bins()),
			count_(parent.bins(), 0),
			parent_(parent)
		{}
		void new_query(unsigned query, Loc seed_offset) {
			last_bin_ = parent_.bin(query / parent_.query_contexts_);
			if (config.trace_pt_membuf) {
				seed_offset_ = seed_offset;
				query_ = query;
				if (buffer_[last_bin_].size() >= buffer_size)
					flush(last_bin_);
			}
			else
				if (text_buffer_[last_bin_].size() >= buffer_size)
					flush(last_bin_);
			
			if (!config.trace_pt_membuf) {
				text_buffer_[last_bin_].write((uint16_t)0);
				text_buffer_[last_bin_].write(query);
				text_buffer_[last_bin_].write(seed_offset);
			}
		}
		void write(unsigned query, PackedLoc subject, uint16_t score, uint32_t target_block_id = 0)
		{
			assert(score > 0);
			++count_[last_bin_];
			assert(last_bin_ < parent_.bins());
			if (config.trace_pt_membuf) {
#ifdef HIT_KEEP_TARGET_ID
				buffer_[last_bin_].emplace_back(query, subject, seed_offset_, score, target_block_id);
#else
				buffer_[last_bin_].emplace_back(query, subject, seed_offset_, score);
#endif
			}
			else {
				text_buffer_[last_bin_].write((uint16_t)score);
				if (parent_.long_subject_offsets_)
					text_buffer_[last_bin_].write_raw((const char*)&subject, 5);
				else
					text_buffer_[last_bin_].write(subject.low);
#ifdef HIT_KEEP_TARGET_ID
				text_buffer_[last_bin_].write(target_block_id);
#endif
			}
		}
		void flush(int bin)
		{
			if(config.trace_pt_membuf) {
				if (buffer_[bin].size() == 0)
					return;
				std::lock_guard<std::mutex> lock(parent_.bin_mutex_[bin]);
				parent_.hit_buf_[bin].insert(parent_.hit_buf_[bin].end(), buffer_[bin].begin(), buffer_[bin].end());
				buffer_[bin].clear();
			}
			else {
				if (text_buffer_[bin].size() == 0)
					return;
				MMap& tmp_file = parent_.tmp_file_[bin];
				const uint32_t size = text_buffer_[bin].size();
				{
					std::lock_guard<std::mutex> lock(parent_.bin_mutex_[bin]);
					tmp_file.write(&size, sizeof(uint32_t));
					tmp_file.write(text_buffer_[bin].data(), text_buffer_[bin].size());
				}
				text_buffer_[bin].clear();
			}			
		}
		virtual ~Writer()
		{
			for (int bin = 0; bin < parent_.bins(); ++bin) {
				flush(bin);
				parent_.count_[bin] += count_[bin];
			}
		}
	private:
		const size_t buffer_size;
		int last_bin_;
		std::vector<std::vector<Hit>> buffer_;
		Loc seed_offset_;
		BlockId query_;
		std::vector<TextBuffer> text_buffer_;
		std::vector<size_t> count_;
		HitBuffer &parent_;
	};

	void load(int64_t max_size);

	std::tuple<Hit*, int64_t, Key, Key> retrieve() {
		if (config.trace_pt_membuf) {
			if (bins_processed_ >= bins())
				return std::tuple<Hit*, int64_t, Key, Key> { nullptr, 0, 0, 0 };
			auto& buf = hit_buf_[bins_processed_];
			++bins_processed_;
			return std::tuple<Hit*, int64_t, Key, Key> { buf.data(), buf.size(), input_range_next_.first, input_range_next_.second };
		}
		else {
			if (load_worker_) {
				load_worker_->join();
				delete load_worker_;
				load_worker_ = nullptr;
			}			
			return std::tuple<Hit*, int64_t, Key, Key> { data_next_, data_size_next_, input_range_next_.first, input_range_next_.second };
		}
	}

	size_t total_disk_size() {
		return total_disk_size_;
	}

	int64_t bin_size(int i) const {
		return count_[i];
	}

	void alloc_buffer() {
		if (config.trace_pt_membuf)
			return;
		int64_t max_size = 0;
		for (int i = 0; i < bins(); ++i)
			max_size = std::max(max_size, bin_size(i));
		alloc_size_ = max_size;
		if (max_size == 0) {
			data_next_ = nullptr;
			return;
		}
#ifdef _MSC_VER
		data_next_ = new Hit[max_size];
#else
		data_next_ = (Hit*)mmap(nullptr, max_size * sizeof(Hit), PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
		if(data_next_ == MAP_FAILED) {
			throw std::runtime_error("Failed to allocate memory for hit buffer");
		}
#endif
	}

	void free_buffer() {
		if (!config.trace_pt_membuf) {
#ifdef _MSC_VER
			delete[] data_next_;
#else
			munmap(data_next_, alloc_size_ * sizeof(Search::Hit));
#endif
		}
	}

private:

	void load_bin(Hit* out, int bin);

	const std::vector<Key> key_partition_;
	const bool long_subject_offsets_;
	const int query_contexts_;
	int bins_processed_;
	int64_t total_disk_size_;
	std::vector<std::vector<Hit>> hit_buf_;
	std::vector<MMap> tmp_file_;
	std::vector<std::mutex> bin_mutex_;
	std::atomic_size_t *count_;
	std::pair<Key, Key> input_range_next_;
	Hit* data_next_;
	int64_t data_size_next_, alloc_size_;
	std::thread* load_worker_;

};

}
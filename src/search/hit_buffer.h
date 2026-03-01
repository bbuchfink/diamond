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
#include <assert.h>
#include <vector>
#include <thread>
#include <atomic>
#include "basic/config.h"
#include "hit.h"
#include "util/io/file.h"
#include "util/data_structures/queue.h"
#include "util/text_buffer.h"
#include "run/config.h"

namespace Search {

struct HitBuffer
{

	using Vector = std::vector<Hit>;
	using Key = uint32_t;
	
	HitBuffer(const std::vector<Key>& key_partition, const std::string& tmpdir, bool long_subject_offsets, int query_contexts, int thread_count, Config& cfg);
	HitBuffer(const std::vector<Key>& key_partition, const std::string& tmpdir, bool long_subject_offsets, int query_contexts, int thread_count, uint32_t max_query, uint64_t max_target, SimpleThreadPool& search_pool);
	~HitBuffer() noexcept(false);
	void finish_writing();

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
			buffer_size(config.trace_pt_membuf ? 8192 : 65535),
			last_bin_(0),
			buffer_(config.trace_pt_membuf ? parent.bins() : 0),
			text_buffer_(config.trace_pt_membuf ? 0 : parent.bins()),
			count_(parent.bins(), 0),
			buf_count_(parent.bins(), 0),
			parent_(parent)
		{
			const int bins = parent.bins();
			if (!config.trace_pt_membuf)
				for (int i = 0; i < bins; ++i)
					text_buffer_[i] = new TextBuffer();
			else
				for (int i = 0; i < bins; ++i) {
					buffer_[i] = new std::vector<Hit>();
					buffer_[i]->reserve(buffer_size);
				}
		}
		void new_query(unsigned query, Loc seed_offset) {
			last_bin_ = parent_.bin(query / parent_.query_contexts_);
			seed_offset_ = seed_offset;
			query_ = query;
			if (!config.trace_pt_membuf) {
				text_buffer_[last_bin_]->write((uint16_t)0);
				text_buffer_[last_bin_]->write(query);
				text_buffer_[last_bin_]->write(seed_offset);
			}
		}
		void write(unsigned query, PackedLoc subject, uint16_t score, uint32_t target_block_id = 0)
		{
			assert(score > 0);			
			assert(last_bin_ < parent_.bins());
			if (config.trace_pt_membuf) {
				if (buffer_[last_bin_]->size() >= buffer_size)
					flush(last_bin_, false);
			}
			else
				if (text_buffer_[last_bin_]->size() >= buffer_size) {
					flush(last_bin_, false);
					text_buffer_[last_bin_]->write((uint16_t)0);
					text_buffer_[last_bin_]->write(query_);
					text_buffer_[last_bin_]->write(seed_offset_);
				}
			++count_[last_bin_];
			++buf_count_[last_bin_];
			if (config.trace_pt_membuf) {
#ifdef HIT_KEEP_TARGET_ID
				buffer_[last_bin_]->emplace_back(query, subject, seed_offset_, score, target_block_id);
#else
				buffer_[last_bin_]->emplace_back(query, subject, seed_offset_, score);
#endif
			}
			else {
				text_buffer_[last_bin_]->write((uint16_t)score);
				if (parent_.long_subject_offsets_)
					text_buffer_[last_bin_]->write_raw((const char*)&subject, 5);
				else
					text_buffer_[last_bin_]->write(subject.low);
#ifdef HIT_KEEP_TARGET_ID
				text_buffer_[last_bin_]->write(target_block_id);
#endif
			}
		}
		void flush(int bin, bool done)
		{
			if(config.trace_pt_membuf) {
				if (buffer_[bin]->size() == 0) {
					delete buffer_[bin];
					return;
				}
				parent_.membuf_out_queue_[bin]->enqueue(std::pair<int, std::vector<Hit>*>(bin, buffer_[bin]));
				if (!done) {
					buffer_[bin] = new std::vector<Hit>();
					buffer_[bin]->reserve(buffer_size);
				}
			}
			else {
				if (text_buffer_[bin]->size() == 0) {
					delete text_buffer_[bin];
					return;
				}
				text_buffer_[bin]->write((uint16_t)0);
				parent_.out_queue_[bin]->enqueue(std::tuple<int, TextBuffer*, uint32_t>(bin, text_buffer_[bin], buf_count_[bin]));
				buf_count_[bin] = 0;
				if(!done) text_buffer_[bin] = new TextBuffer();
			}
		}
		virtual ~Writer()
		{
			for (int bin = 0; bin < parent_.bins(); ++bin) {
				flush(bin, true);
				parent_.count_[bin] += count_[bin];
			}
		}
	private:
		const size_t buffer_size;
		int last_bin_;
		std::vector<std::vector<Hit>*> buffer_;
		Loc seed_offset_;
		BlockId query_;
		std::vector<TextBuffer*> text_buffer_;
		std::vector<size_t> count_;
		std::vector<uint32_t> buf_count_;
		HitBuffer &parent_;
	};

	bool load(size_t max_size);

	std::tuple<Hit*, size_t, Key, Key> retrieve() {
		if (config.trace_pt_membuf || config.swipe_all) {
			if (bins_processed_ >= bins())
				return std::tuple<Hit*, size_t, Key, Key> { nullptr, 0, 0, 0 };
			++bins_processed_;
			if(config.swipe_all)
				return std::tuple<Hit*, size_t, Key, Key> { nullptr, 0, input_range_next_.first, input_range_next_.second };
			else {
				auto& buf = hit_buf_[bins_processed_];
				return std::tuple<Hit*, size_t, Key, Key> { buf.data(), buf.size(), input_range_next_.first, input_range_next_.second };
			}
		}
		else {
			if (load_worker_) {
				load_worker_->join();
				delete load_worker_;
				load_worker_ = nullptr;
				if(load_exception_)
					std::rethrow_exception(load_exception_);
			} else
				throw std::runtime_error("HitBuffer retrieve w/o load");
			std::swap(data_loading_, data_finished_);
			std::swap(mmap_finished_, mmap_loading_);
			return std::tuple<Hit*, size_t, Key, Key> { data_finished_, data_size_next_, input_range_next_.first, input_range_next_.second };
		}
	}

	size_t total_disk_size() {
		return total_disk_size_;
	}

	int64_t bin_size(int i) const {
		return count_[i];
	}

	uint64_t next_bin_size() const {
		return bins_processed_ < bins() ? bin_size(bins_processed_) : 0;
	}

	void alloc_buffer();
	void free_buffer();

private:

	void load_bin(Hit* out, int bin);
	void write_worker(const std::atomic<bool>& stop, int bin);

	const std::vector<Key> key_partition_;
	const bool long_subject_offsets_;
	const int query_contexts_;
	const uint32_t max_query_;
	const uint64_t max_target_;
	SimpleThreadPool& search_pool_;
	int bins_processed_;
	int64_t total_disk_size_;
	std::vector<std::vector<Hit>> hit_buf_;
	std::vector<File> tmp_file_;
	std::atomic_size_t *count_;
	std::pair<Key, Key> input_range_next_;
	Hit* data_loading_, *data_finished_;
	bool mmap_loading_, mmap_finished_;
	int64_t data_size_next_, alloc_size_;
	std::thread* load_worker_;
	std::exception_ptr load_exception_;
	std::vector<std::thread::id> writer_;
	std::vector<Queue<std::tuple<int, TextBuffer*, uint32_t>>*> out_queue_;
	std::vector<Queue<std::pair<int, std::vector<Hit>*>>*> membuf_out_queue_;

};

}
/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

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
#include <assert.h>
#include <vector>
#include <thread>
#include <atomic>
#include <mutex>
#include "basic/config.h"
#include "hit.h"
#include "util/io/temp_file.h"
#include "util/io/compressed_buffer.h"

namespace Search {

struct HitBuffer
{

	using Vector = std::vector<Hit>;
	using Key = uint32_t;
	
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
			seed_offset_ = seed_offset;
			query_ = query;
			if (config.trace_pt_membuf) {				
				if (buffer_[last_bin_].size() >= buffer_size)
					flush(last_bin_);
			}
			else
				if (text_buffer_[last_bin_].size() >= buffer_size)
					flush(last_bin_);
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
				Hit h;
				h.query_ = query;
				h.subject_ = subject;
				h.seed_offset_ = seed_offset_;
				h.score_ = score;				
#ifdef HIT_KEEP_TARGET_ID
				h.target_block_id = target_block_id;
#endif
				text_buffer_[last_bin_].write(h);
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
				text_buffer_[bin].finish();
				if (text_buffer_[bin].size() == 0)
					return;
				TmpFile& tmp_file = parent_.tmp_file_[bin];
				{
					std::lock_guard<std::mutex> lock(parent_.bin_mutex_[bin]);
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
		std::vector<CompressedBuffer> text_buffer_;
		std::vector<size_t> count_;
		HitBuffer &parent_;
	};

	void load(size_t max_size);

	std::tuple<Hit*, size_t, Key, Key> retrieve() {
		if (config.trace_pt_membuf) {
			if (bins_processed_ >= bins())
				return std::tuple<Hit*, size_t, Key, Key> { nullptr, 0, 0, 0 };
			auto& buf = hit_buf_[bins_processed_];
			++bins_processed_;
			return std::tuple<Hit*, size_t, Key, Key> { buf.data(), buf.size(), input_range_next_.first, input_range_next_.second };
		}
		else {
			if (load_worker_) {
				load_worker_->join();
				delete load_worker_;
				load_worker_ = nullptr;
			}			
			return std::tuple<Hit*, size_t, Key, Key> { data_next_, data_size_next_, input_range_next_.first, input_range_next_.second };
		}
	}

	size_t total_disk_size() {
		return total_disk_size_;
	}

	int64_t bin_size(int i) const {
		return count_[i];
	}

	void alloc_buffer();
	void free_buffer();

private:

	void load_bin(Hit* out, int bin);

	const std::vector<Key> key_partition_;
	const bool long_subject_offsets_;
	const int query_contexts_;
	int bins_processed_;
	int64_t total_disk_size_;
	std::vector<std::vector<Hit>> hit_buf_;
	std::vector<TmpFile> tmp_file_;
	std::vector<std::mutex> bin_mutex_;
	std::atomic_size_t *count_;
	std::pair<Key, Key> input_range_next_;
	Hit* data_next_;
	bool mmap_;
	int64_t data_size_next_, alloc_size_;
	std::thread* load_worker_;

};

}
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

#ifndef WIN32
#include <sys/mman.h>
#include <fcntl.h>
#endif
#include "hit_buffer.h"
#include "basic/config.h"
#include "util/log_stream.h"
#include "util/io/output_file.h"
#include "util/parallel/simple_thread_pool.h"

using std::vector;
using std::string;
using std::thread;
using std::copy;
using std::runtime_error;
using std::mutex;
using std::endl;
using std::atomic_size_t;
using std::pair;
using std::atomic;
using std::tuple;

namespace Search {

HitBuffer::HitBuffer(const vector<Key>& key_partition, const string& tmpdir, bool long_subject_offsets, int query_contexts, int thread_count) :
	key_partition_(key_partition),
	long_subject_offsets_(long_subject_offsets),
	query_contexts_(query_contexts),
	bins_processed_(0),
	total_disk_size_(0),
	mmap_(false),
	load_worker_(nullptr)
{
	log_stream << "Async_buffer() " << key_partition.back() << std::endl;
	count_ = new atomic_size_t[key_partition.size()];
	for (size_t i = 0; i < key_partition.size(); ++i) {
		if (config.trace_pt_membuf) {
			membuf_out_queue_.push_back(new Queue<std::pair<int, std::vector<Hit>*>>(thread_count * 4, thread_count, 1, pair<int, vector<Hit>*>(0, nullptr)));			
			hit_buf_.emplace_back();
			hit_buf_.back().reserve(25 * GIGABYTES / sizeof(Hit));			
		}
		else {
			out_queue_.push_back(new Queue<tuple<int, TextBuffer*, uint32_t>>(thread_count * 4, thread_count, 1, tuple<int, TextBuffer*, uint32_t>(0, nullptr, 0)));
			tmp_file_.push_back(File(Temporary()));
		}
		writer_.push_back(new thread(&HitBuffer::write_worker, this, (int)i));
		count_[i].store(0, std::memory_order_relaxed);
	}	
}

void HitBuffer::write_worker(int bin) {
	try {
		if (config.trace_pt_membuf) {
			for (;;) {
				pair<int, vector<Hit>*> buf;
				if (!membuf_out_queue_[bin]->wait_and_dequeue(buf))
					break;
				vector<Hit>& hits = hit_buf_[buf.first];
				hits.insert(hits.end(), buf.second->begin(), buf.second->end());
				delete buf.second;
			}
		}
		else {
			for (;;) {
				tuple<int, TextBuffer*, uint32_t> buf;
				if (!out_queue_[bin]->wait_and_dequeue(buf))
					break;
				File& tmp_file = tmp_file_[std::get<0>(buf)];
				//std::lock_guard<std::mutex> lock(bin_mutex_[bin]);
				tmp_file.write(std::get<1>(buf)->size());
				tmp_file.write(std::get<2>(buf));
				tmp_file.write(std::get<1>(buf)->data(), std::get<1>(buf)->size());
				delete std::get<1>(buf);
			}
		}
	}
	catch (...) {

	}
}

void HitBuffer::finish_writing() {
	for (int i = 0; i < out_queue_[0]->producer_count(); ++i)
		if (config.trace_pt_membuf) {
			for (int i = 0; i < bins(); ++i)
				membuf_out_queue_[i]->close();
		}
		else {
			for (int i = 0; i < bins(); ++i)
				out_queue_[i]->close();
		}
	for (auto& t : writer_) {
		t->join();
		delete t;
	}
	for(auto& q : membuf_out_queue_)
		delete q;
	for (auto& q : out_queue_)
		delete q;
	writer_.clear();
	membuf_out_queue_.clear();
	out_queue_.clear();
}

HitBuffer::~HitBuffer() {
	if (!writer_.empty())
		throw runtime_error("HitBuffer::~HitBuffer(): writer thread still active");
	delete[] count_;
}

bool HitBuffer::load(size_t max_size) {
	max_size = std::max(max_size, (size_t)1);
	data_size_next_ = 0;
	auto worker = [this](int end) {
		Hit* out = data_next_;
		for (; bins_processed_ < end; ++bins_processed_) {
			load_bin(out, bins_processed_);
			out += count_[bins_processed_];
		}
		};
	if (bins_processed_ == bins()) {
		return false;
	}
	size_t size = count_[bins_processed_], current_size;
	const int begin = bins_processed_;
	int end = bins_processed_ + 1;
	if (!config.trace_pt_membuf) {
		size_t disk_size = tmp_file_[bins_processed_].size();
		// consider using more bins here
		while (end < bins() && (size + (current_size = count_[end])) * sizeof(Hit) < max_size && (end - bins_processed_ == 0)) {
			size += current_size;
			disk_size += tmp_file_[end].size();
			++end;
		}
		log_stream << "Async_buffer.load() " << size << " (" << (double)size * sizeof(Hit) / (1 << 30) << " GB, " << (double)disk_size / (1 << 30) << " GB on disk)" << endl;
		total_disk_size_ += disk_size;
		data_size_next_ = size;
		load_worker_ = new thread(worker, end);
	}
	input_range_next_.first = this->begin(begin);
	input_range_next_.second = this->end(end - 1);	
	return true;
}

void HitBuffer::load_bin(Hit* out, int bin)
{
	if (config.trace_pt_membuf)
		return;
#if !_MSC_VER && !__APPLE__
	if (bin < bins() - 1)
		posix_fadvise(fileno(tmp_file_[bin + 1].file()), 0, 0, POSIX_FADV_SEQUENTIAL | POSIX_FADV_WILLNEED);
#endif
	File& f = tmp_file_[bin];
	if (count_[bin] == 0) {
		f.close();
		return;
	}
	const int p = config.threads_ > 1 ? std::min(config.threads_ - 1, config.load_threads) : 1;
	Queue<pair<vector<char>*, uint32_t>> queue(p * 4, 1, p, pair<vector<char>*, uint32_t>(nullptr, 0));
	atomic<Hit*> out_ptr = out;
	atomic<uint64_t> count = 0;
	f.seek(0, SEEK_SET);
	SimpleThreadPool pool;
	auto worker = [&](const atomic<bool>& stop) {
		uint64_t my_count = 0;
		while (!stop) {
			pair<vector<char>*, uint32_t> v;
			if (!queue.wait_and_dequeue(v))
				break;
			Hit* dst = out_ptr.fetch_add(v.second, std::memory_order_relaxed);
			vector<char>::const_iterator ptr = v.first->begin(), end = v.first->end();
			uint16_t nullscore;
			memcpy(&nullscore, &*ptr, 2);
			ptr += 2;
			while(ptr < end) {
				uint32_t query_id, seed_offset;
				memcpy(&query_id, &*ptr, 4);
				ptr += 4;
				memcpy(&seed_offset, &*ptr, 4);
				ptr += 4;
				PackedLoc subject_loc;
				uint32_t x;
				while(ptr<end) {
					uint16_t score;
					memcpy(&score, &*ptr, 2);
					ptr += 2;
					if (score == 0)
						break;
					if (long_subject_offsets_) {
						memcpy(&subject_loc, &*ptr, sizeof(PackedLoc));
						ptr += sizeof(PackedLoc);
					}
					else {
						memcpy(&x, &*ptr, 4);
						ptr += 4;
						subject_loc = x;
					}
#ifdef HIT_KEEP_TARGET_ID
					uint32_t target_block_id;
					memcpy(&target_block_id, &*ptr, 4);
					ptr += 4;
					dst->target_block_id = target_block_id;
#endif
					dst->query_ = query_id;
					dst->subject_ = subject_loc;
					dst->seed_offset_ = seed_offset;
					dst->score_ = score;
					++dst;
					++my_count;
				}
			}
		}
		count.fetch_add(my_count, std::memory_order_relaxed);
		};
	for (int i = 0; i < p; ++i)
		pool.spawn(worker);
	while(!pool.stop()) {
		uint32_t n;
		size_t l;
		if (f.read_max(&l, sizeof(size_t)) < sizeof(size_t))
			break;
		f.read(n);
		vector<char>* v = new vector<char>(l);
		f.read(v->data(), l);
		queue.enqueue(pair<vector<char>*, uint32_t>(v, n));
	}
	queue.close();
	pool.join_all();

	if (count != count_[bin])
		throw runtime_error("Mismatching hit count / possibly corrupted temporary file");
	f.close();
}

void HitBuffer::alloc_buffer() {
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
	mmap_ = false;
#else
	int flags = MAP_PRIVATE | MAP_ANONYMOUS;
#ifdef MAP_HUGETLB
	flags |= MAP_HUGETLB;
#endif
	data_next_ = (Hit*)mmap(nullptr, max_size * sizeof(Hit), PROT_READ | PROT_WRITE, flags, -1, 0);
	if (data_next_ == MAP_FAILED) {
		data_next_ = new Hit[max_size];
		mmap_ = false;
	}
	else
		mmap_ = true;
#endif
}

void HitBuffer::free_buffer() {
	if (!config.trace_pt_membuf) {
#ifdef _MSC_VER
		delete[] data_next_;
#else
		if (mmap_)
			munmap(data_next_, alloc_size_ * sizeof(Search::Hit));
		else
			delete[] data_next_;
#endif
	}

}

}
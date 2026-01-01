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

using std::vector;
using std::string;
using std::thread;
using std::copy;
using std::runtime_error;
using std::mutex;
using std::endl;
using std::atomic_size_t;

namespace Search {

HitBuffer::HitBuffer(const vector<Key>& key_partition, const string& tmpdir, bool long_subject_offsets, int query_contexts) :
	key_partition_(key_partition),
	long_subject_offsets_(long_subject_offsets),
	query_contexts_(query_contexts),
	bins_processed_(0),
	total_disk_size_(0),
	bin_mutex_(key_partition.size()),
	mmap_(false),
	load_worker_(nullptr)
{
	log_stream << "Async_buffer() " << key_partition.back() << std::endl;
	count_ = new atomic_size_t[key_partition.size()];
	for (size_t i = 0; i < key_partition.size(); ++i) {
		if (config.trace_pt_membuf) {
			hit_buf_.emplace_back();
			hit_buf_.back().reserve(GIGABYTES / sizeof(Hit));
		}
		else {
			tmp_file_.emplace_back();
		}
		count_[i].store(0, std::memory_order_relaxed);
	}
}

bool HitBuffer::load(size_t max_size) {
	max_size = std::max(max_size, (size_t)1);
	auto worker = [this](int end) {
		Hit* out = data_next_;
		for (; bins_processed_ < end; ++bins_processed_) {
			load_bin(out, bins_processed_);
			out += count_[bins_processed_];
		}
		};
	if (bins_processed_ == bins()) {
		data_next_ = nullptr;
		return false;
	}
	size_t size = count_[bins_processed_], current_size;
	const int begin = bins_processed_;
	int end = bins_processed_ + 1;
	if (!config.trace_pt_membuf) {
		size_t disk_size = tmp_file_[bins_processed_].size();
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
	if (!config.trace_pt_membuf) {
#if !_MSC_VER && !__APPLE__
		if (bin < bins() - 1)
			posix_fadvise(fileno(tmp_file_[bin + 1].file()), 0, 0, POSIX_FADV_SEQUENTIAL | POSIX_FADV_WILLNEED);
#endif
#ifdef WITH_ZSTD
		const Compressor c = Compressor::ZSTD;
#else
		const Compressor c = Compressor::ZLIB;
#endif
		if (count_[bin] > 0) {
			tmp_file_[bin].seek(0, SEEK_SET);
			size_t n = decompress(tmp_file_[bin].file(), out, count_[bin] * sizeof(Hit), c);
			if (n != count_[bin] * sizeof(Hit))
				throw runtime_error("Mismatching hit count / possibly corrupted temporary file");
		}
		tmp_file_[bin].close();
	}
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
#else
	int flags = MAP_PRIVATE | MAP_ANONYMOUS;
#ifdef MAP_HUGETLB
	flags |= MAP_HUGETLB;
#endif
	data_next_ = (Hit*)mmap(nullptr, max_size * sizeof(Hit), PROT_READ | PROT_WRITE, flags, -1, 0);
	if (data_next_ == MAP_FAILED) {
		data_next_ = new Hit[max_size];
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
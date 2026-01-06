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
#include "external.h"
#include "util/io/input_file.h"
#include "util/log_stream.h"
#include "util/parallel/atomic.h"
#include "util/string/string.h"

namespace Cluster {

template<typename T>
std::vector<std::string> radix_cluster(Job& job, const VolumedFile& bucket, const std::string& output_dir, int bits_unsorted) {
	const int64_t BUF_SIZE = 4096;
	const uint64_t shift = bits_unsorted - RADIX_BITS;

	std::unique_ptr<FileArray> output_files(new FileArray(output_dir, RADIX_COUNT, job.worker_id()));
	std::atomic<int64_t> next(0);

	std::vector<std::thread> workers;
	auto worker = [&](int thread_id) {
		BufferArray buffers(*output_files, RADIX_COUNT);
		int64_t v = 0;
		while (v = next.fetch_add(1, std::memory_order_relaxed), v < (int64_t)bucket.size()) {
			InputFile in(bucket[v].path);
			const int64_t n = bucket[v].record_count;
			std::unique_ptr<T[]> data(new T[n]);
			in.read(data.get(), n);
			in.close();
			const T* end = data.get() + n;
			for (const T* ptr = data.get(); ptr < end; ++ptr) {
				const int radix = (ptr->key() >> shift) & (RADIX_COUNT - 1);
				buffers.write(radix, *ptr);
			}
		}
		};
	for (int i = 0; i < std::min(config.threads_, (int)bucket.size()); ++i)
		workers.emplace_back(worker, i);
	for (auto& t : workers)
		t.join();
	std::vector<std::string> buckets;
	buckets.reserve(RADIX_COUNT);
	for (int i = 0; i < RADIX_COUNT; ++i)
		buckets.push_back(output_files->bucket(i));
	TaskTimer timer("Closing the output files");
	output_files.reset();
	timer.finish();
	job.log("Radix sorted bucket records=%zu", bucket.records());
	return buckets;
}

template<typename T>
std::vector<std::string> radix_sort(Job& job, const std::vector<std::string>& buckets, int bits_unsorted) {
	const std::string base_path = Cluster::base_path(buckets.front()),
		queue_path = base_path + PATH_SEPARATOR + "radix_sort_queue",
		result_path = base_path + PATH_SEPARATOR + "radix_sort_out";
	const uint64_t size_limit = Util::String::interpret_number(config.memory_limit.get(DEFAULT_MEMORY_LIMIT));
	Atomic queue(queue_path);
	FileStack out(result_path);
	int64_t i, buckets_processed = 0;
	while (i = queue.fetch_add(), i < (int64_t)buckets.size()) {
		VolumedFile bucket(buckets[i]);
		const size_t data_size = bucket.records() * sizeof(T);
		job.log("Radix sorting. Bucket=%lli/%lli Records=%s Size=%s", i + 1, buckets.size(), Util::String::format(bucket.records()).c_str(), Util::String::format(data_size).c_str());
		if (data_size > size_limit) {
			const std::vector<std::string> v = radix_cluster<T>(job, bucket, path(buckets[i]), bits_unsorted);
			out.lock();
			for (const std::string& s : v)
				out.push_non_locked(s);
			out.unlock();
		}
		else if (bucket.records() > 0)
			out.push(buckets[i]);
		else
			bucket.remove();
		++buckets_processed;
	}
	Atomic finished(base_path + PATH_SEPARATOR + "radix_sort_finished");
	finished.fetch_add(buckets_processed);
	finished.await(buckets.size());
	return read_list(result_path);
}

}
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
#include <atomic>
#include "external.h"
#include "util/io/input_file.h"
#include "util/log_stream.h"
#include "util/parallel/atomic.h"
#include "util/string/string.h"

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
	job.log("Radix sorted bucket records=%zu", bucket.sparse_records());
	return buckets;
}

template<typename T>
std::vector<std::string> radix_sort(Job& job, const std::vector<std::string>& buckets, int bits_unsorted) {
	const std::string base_path = ::base_path(buckets.front()),
		queue_path = base_path + PATH_SEPARATOR + "radix_sort_queue",
		result_path = base_path + PATH_SEPARATOR + "radix_sort_out";
	const uint64_t size_limit = Util::String::interpret_number(config.memory_limit.get(DEFAULT_MEMORY_LIMIT));
	Atomic queue(queue_path);
	FileStack out(result_path);
	int64_t i, buckets_processed = 0;
	while (i = queue.fetch_add(), i < (int64_t)buckets.size()) {
		VolumedFile bucket(buckets[i]);
		const size_t data_size = bucket.sparse_records() * sizeof(T);
		job.log("Radix sorting. Bucket=%lli/%lli Records=%s Size=%s", i + 1, buckets.size(), Util::String::format(bucket.sparse_records()).c_str(), Util::String::format(data_size).c_str());
		if (data_size > size_limit) {
			const std::vector<std::string> v = radix_cluster<T>(job, bucket, path(buckets[i]), bits_unsorted);
			out.lock();
			for (const std::string& s : v)
				out.push_non_locked(s);
			out.unlock();
		}
		else if (bucket.sparse_records() > 0)
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
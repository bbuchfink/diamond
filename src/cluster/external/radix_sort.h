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
#include <type_traits>
#include <sstream>
#include "external.h"
#include "util/io/input_file.h"
#include "util/log_stream.h"
#include "util/parallel/atomic.h"
#include "util/parallel/simple_thread_pool.h"
#include "util/string/string.h"
#include "file_array.h"
#include "util/memory/memory_resource.h"

template<typename T>
typename std::enable_if<T::POD>::type
radix_cluster_read(InputFile& in, size_t n, BufferArray& buffers, uint64_t shift) {
	T x;
	for (size_t i = 0; i < n; ++i) {
 		in.read(&x, sizeof(T));
		const uint64_t radix = (x.key() >> shift) & (RADIX_COUNT - 1);
		buffers.write(radix, x);
	}
}

template<typename T>
typename std::enable_if<!T::POD>::type
radix_cluster_read(InputFile& in, size_t n, BufferArray& buffers, uint64_t shift) {
	std::pmr::monotonic_buffer_resource pool;
	T x(pool);
	for (size_t i = 0; i < n; ++i) {
 		deserialize(in, x);
		const uint64_t radix = (x.key() >> shift) & (RADIX_COUNT - 1);
		buffers.write(radix, x);
	}
}

template<typename T>
RadixedTable radix_cluster(Job& job, const VolumedFile& bucket, const std::string& output_dir, int bits_unsorted) {
	const int64_t BUF_SIZE = 4096;
	const uint64_t shift = bits_unsorted - RADIX_BITS;

	std::unique_ptr<FileArray> output_files(new FileArray(output_dir, RADIX_COUNT, job.worker_id(), true));
	std::atomic<int64_t> next(0);

	SimpleThreadPool pool;
	auto worker = [&](const std::atomic<bool>& stop, int thread_id) {
		BufferArray buffers(*output_files, RADIX_COUNT);
		int64_t v = 0;

		while (!stop.load(std::memory_order_relaxed) && (v = next.fetch_add(1, std::memory_order_relaxed), v < (int64_t)bucket.size())) {
			InputFile in(bucket[v].path);
			const size_t n = bucket[v].record_count;
			radix_cluster_read<T>(in, n, buffers, shift);
			in.close();
		}
		};
	for (int i = 0; i < std::min(config.threads_, (int)bucket.size()); ++i)
		pool.spawn(worker, i);
	pool.join_all();
	TaskTimer timer("Closing the output files");
	const RadixedTable out = output_files->buckets();
	output_files.reset();
	timer.finish();
	job.log("Radix sorted bucket records=%zu", bucket.sparse_records());
	return out;
}

template<typename T>
RadixedTable radix_sort(Job& job, const RadixedTable& buckets, int bits_unsorted) {
	if (bits_unsorted < RADIX_BITS)
		return buckets;
	const std::string base_path = ::base_path(buckets.front().path),
		queue_path = base_path + PATH_SEPARATOR + "radix_sort_queue",
		result_path = base_path + PATH_SEPARATOR + "radix_sort_out";
	Atomic queue(queue_path);
	FileStack out(result_path);
	int64_t i, buckets_processed = 0;
	while (i = queue.fetch_add(), i < (int64_t)buckets.size()) {
		VolumedFile bucket(buckets[i]);
		const size_t data_size = bucket.sparse_records() * sizeof(T);
		job.log("Radix sorting. Bucket=%lli/%lli Records=%s Size=%s", i + 1, buckets.size(), Util::String::format(bucket.sparse_records()).c_str(), Util::String::format(data_size).c_str());
		if (data_size > job.mem_limit) {
			const RadixedTable v = radix_cluster<T>(job, bucket, buckets[i].containing_directory(), bits_unsorted);
			v.append(out);
		}
		else if (bucket.sparse_records() > 0) {
			std::ostringstream ss;
			ss << buckets[i].path << '\t' << bucket.sparse_records() << std::endl;
			out.push(ss.str());
		}
		else
			bucket.remove();
		++buckets_processed;
	}
	Atomic finished(base_path + PATH_SEPARATOR + "radix_sort_finished");
	finished.fetch_add(buckets_processed);
	finished.await(buckets.size());
	return RadixedTable(result_path);
}

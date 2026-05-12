/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

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

namespace External {

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
	const std::string base_path = ::External::base_path(buckets.front().path),
		queue_path = base_path + PATH_SEPARATOR + "radix_sort_queue",
		result_path = base_path + PATH_SEPARATOR + "radix_sort_out";
	Atomic queue(queue_path, job);
	FileStack out(result_path, job);
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
	Atomic finished(base_path + PATH_SEPARATOR + "radix_sort_finished", job);
	finished.fetch_add(buckets_processed);
	finished.await(buckets.size());
	return RadixedTable(result_path);
}

}
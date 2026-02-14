/****
Copyright � 2012-2026 Benjamin J. Buchfink <buchfink@gmail.com>

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

#ifdef _WIN32
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif
#define _REENTRANT
#include "lib/ips4o/ips4o.hpp"
#include "build_pair_table.h"
#include "external.h"
#include "util/parallel/simple_thread_pool.h"

using std::atomic;
using std::vector;
using std::string;

RadixedTable build_pair_table(Job& job, const RadixedTable& seed_table, int shape, int64_t max_oid, FileArray& output_files) {
	const int64_t BUF_SIZE = 4096, shift = bit_length(max_oid) - RADIX_BITS; // promiscuous_cutoff = db_size / config.promiscuous_seed_ratio;
	const string seed_table_base = job.base_dir() + PATH_SEPARATOR + "seed_table_" + std::to_string(shape);
	const string queue_path = seed_table_base + PATH_SEPARATOR + "build_pair_table_queue";
	const bool unid = !config.mutual_cover.present();
	Atomic queue(queue_path);
	atomic<size_t> buckets_processed(0);
	SimpleThreadPool pool;
	const int concurrent_buckets = std::min((int)seed_table.max_buckets(job.mem_limit, sizeof(SeedEntry)), config.threads_);
	const int bucket_workers = div_up(config.threads_, (int)concurrent_buckets);
	job.log("Building pair table. Concurrent buckets=%i Workers per bucket=%i", concurrent_buckets, bucket_workers);
	auto worker = [&](const atomic<bool>& stop) {
		int64_t bucket;
		while (!stop.load(std::memory_order_relaxed) && (bucket = queue.fetch_add(), bucket < (int64_t)seed_table.size())) {
			VolumedFile file(seed_table[bucket]);
			InputBuffer<SeedEntry> data(file, bucket_workers);
			job.log("Building pair table. Bucket=%lli/%lli Records=%s Size=%s", bucket + 1, seed_table.size(), Util::String::format(data.size()).c_str(), Util::String::format(data.byte_size()).c_str());
			ips4o::parallel::sort(data.begin(), data.end(), std::less<SeedEntry>(), bucket_workers);
			auto bucket_worker = [&](const atomic<bool>& stop, int thread_id) {
				BufferArray buffers(output_files, RADIX_COUNT);
				auto it = merge_keys(data.begin(thread_id), data.end(thread_id), SeedEntry::Key());
				while (!stop.load(std::memory_order_relaxed) && it.good()) {
					/*if (it.count() >= promiscuous_cutoff) {
						++it;
						continue;
					}*/
					if (unid)
						get_pairs_uni_cov(it, buffers);
					else
						get_pairs_mutual_cov(it, buffers);
					++it;
				}
				};
			vector<std::thread::id> threads;
			for (int i = 0; i < data.parts(); ++i)
				threads.push_back(pool.spawn(bucket_worker, i));
			pool.join(threads.begin(), threads.end());
			file.remove();
			buckets_processed.fetch_add(1, std::memory_order_relaxed);
		}
		};
	vector<std::thread::id> threads;
	for (int i = 0; i < concurrent_buckets; ++i)
		threads.push_back(pool.spawn(worker));
	pool.join(threads.begin(), threads.end());
	const RadixedTable buckets = output_files.buckets();
	Atomic finished(seed_table_base + PATH_SEPARATOR + "pair_table_finished");
	finished.fetch_add(buckets_processed);
	finished.await(seed_table.size());
	return buckets;
}
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

#define _REENTRANT
#include "lib/ips4o/ips4o.hpp"
#include "build_pair_table.h"
#include "external.h"
#include "util/parallel/simple_thread_pool.h"

using std::atomic;
using std::vector;
using std::string;

namespace External {

RadixedTable build_pair_table(Job& job, const RadixedTable& seed_table, int shape, int64_t max_oid, FileArray& output_files) {
	const int64_t BUF_SIZE = 4096, shift = bit_length(max_oid) - RADIX_BITS; // promiscuous_cutoff = db_size / config.promiscuous_seed_ratio;
	const string seed_table_base = job.base_dir() + PATH_SEPARATOR + "seed_table_" + std::to_string(shape);
	const string queue_path = seed_table_base + PATH_SEPARATOR + "build_pair_table_queue";
	const bool unid = !config.mutual_cover.present();
	Atomic queue(queue_path, job);
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
	Atomic finished(seed_table_base + PATH_SEPARATOR + "pair_table_finished", job);
	finished.fetch_add(buckets_processed);
	finished.await(seed_table.size());
	return buckets;
}

}
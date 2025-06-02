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
		std::array<std::vector<T>, RADIX_COUNT> buffers;
		for (int64_t i = 0; i < RADIX_COUNT; ++i)
			buffers[i].reserve(BUF_SIZE);
		int64_t v = 0;
		while (v = next.fetch_add(1, std::memory_order_relaxed), v < (int64_t)bucket.size()) {
			InputFile in(bucket[v].path, InputFile::NO_AUTODETECT);
			const int64_t n = bucket[v].record_count;
			std::unique_ptr<T[]> data(new T[n]);
			in.read(data.get(), n);
			in.close();
			const T* end = data.get() + n;
			for (const T* ptr = data.get(); ptr < end; ++ptr) {
				const uint64_t radix = (ptr->key() >> shift) & (RADIX_COUNT - 1);
				std::vector<T>& buf = buffers[radix];
				buf.push_back(*ptr);
				if (buf.size() >= BUF_SIZE) {
					//output_files->write(radix, buf.data(), buf.size());
					buf.clear();
				}
			}
		}
		for (int i = 0; i < RADIX_COUNT; ++i)
			;// output_files->write(i, buffers[i].data(), buffers[i].size());
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
	return buckets;
}

template<typename T>
std::vector<std::string> radix_sort(Job& job, const std::vector<std::string>& buckets, int bits_unsorted) {
	const std::string base_path = Cluster::base_path(buckets.front()),
		queue_path = base_path + PATH_SEPARATOR + "radix_sort_queue",
		result_path = base_path + PATH_SEPARATOR + "radix_sort_out";
	const int64_t size_limit = Util::String::interpret_number(config.memory_limit.get(DEFAULT_MEMORY_LIMIT));
	Atomic queue(queue_path);
	FileStack out(result_path);
	int64_t i, buckets_processed = 0;
	while (i = queue.fetch_add(), i < (int64_t)buckets.size()) {
		VolumedFile bucket(buckets[i]);
		const int64_t data_size = bucket.records() * sizeof(T);
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
		++buckets_processed;
	}
	Atomic finished(base_path + PATH_SEPARATOR + "radix_sort_finished");
	finished.fetch_add(buckets_processed);
	finished.await(buckets.size());
	return read_list(result_path);
}

}
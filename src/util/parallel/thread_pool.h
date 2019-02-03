#ifndef THREAD_POOL_H_
#define THREAD_POOL_H_

#include <thread>
#include <atomic>
#include <vector>

namespace Util { namespace Parallel {

template<typename _f, typename... _args>
void pool_worker(std::atomic<size_t> *partition, size_t thread_id, size_t partition_count, _f f, _args... args) {
	size_t p;
	while ((p = (*partition)++) < partition_count) {
		f(p, thread_id, args...);
	}
}

template<typename _f, typename... _args>
void scheduled_thread_pool(size_t thread_count, _f f, _args... args) {
	std::atomic<size_t> partition(0);
	std::vector<std::thread> threads;
	for (size_t i = 0; i < thread_count; ++i)
		threads.emplace_back(f, &partition, i, args...);
	for (std::thread &t : threads)
		t.join();
}

template<typename _f, typename... _args>
void scheduled_thread_pool_auto(size_t thread_count, size_t partition_count, _f f, _args... args) {
	scheduled_thread_pool(thread_count, pool_worker<_f, _args...>, partition_count, f, args...);
}

}}

#endif
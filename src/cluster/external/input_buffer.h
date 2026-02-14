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
#include "util/memory/memory_resource.h"
#include "volume.h"
#include "util/algo/partition.h"
#include "util/parallel/simple_thread_pool.h"

template<typename T>
struct InputBuffer {

	using Container = std::conditional_t<T::POD, std::vector<T>, std::pmr::list<T>>;
	using Iterator = typename Container::iterator;
	using ConstIterator = typename Container::const_iterator;

	InputBuffer(const VolumedFile& f, int parts = config.threads_) :
		size_(f.sparse_records()),
		data_(size_),
		part_(size_, parts)
	{
		std::atomic<int64_t> next(0);
		SimpleThreadPool pool;
		auto worker = [&](const std::atomic<bool>& stop) {
			int64_t v;
			while (!stop.load(std::memory_order_relaxed) && (v = next.fetch_add(1, std::memory_order_relaxed), v < (int64_t)f.size())) {
				InputFile in(f[v].path);
				in.read(&data_[f[v].oid_begin], f[v].record_count);
				in.close();
			}
			};
		for (int i = 0; i < std::min(config.threads_, (int)f.size()); ++i)
			pool.spawn(worker);
		pool.join_all();
	}

	InputBuffer(const VolumedFile& f, std::pmr::monotonic_buffer_resource& mem_pool, int parts = config.threads_) :
		size_(f.sparse_records()),
		data_(&mem_pool),
		part_(size_, parts)
	{
		std::atomic<int64_t> next(0);
		SimpleThreadPool pool;
		std::mutex mtx;
		auto worker = [&](const std::atomic<bool>& stop) {
			int64_t v;
			while (!stop.load(std::memory_order_relaxed) && (v = next.fetch_add(1, std::memory_order_relaxed), v < (int64_t)f.size())) {
				std::pmr::list<T> lst(&mem_pool);
				InputFile in(f[v].path);
				T x(mem_pool);
				for (size_t i = 0; i < f[v].record_count; ++i) {
					deserialize(in, x);
					lst.push_back(std::move(x));
				}
				in.close();
				{
					std::lock_guard<std::mutex> lock(mtx);
					data_.splice(data_.end(), lst);
				}
			}
			};
		for (int i = 0; i < std::min(config.threads_, (int)f.size()); ++i)
			pool.spawn(worker);
		pool.join_all();
	}

	size_t size() const {
		return size_;
	}

	size_t byte_size() const {
		return size_ * sizeof(T);
	}

	Iterator begin() {
		return data_.begin();
	}

	Iterator end() {
		return data_.end();
	}

	ConstIterator cbegin() const {
		return data_.cbegin();
	}

	ConstIterator cend() const {
		return data_.cend();
	}

	ConstIterator begin(int part) const {
		ConstIterator begin = data_.cbegin() + part_.begin(part), end = this->cend();
		if (part > 0)
			while (begin < end && begin[-1].key() == begin[0].key())
				++begin;
		return begin;
	}

	ConstIterator end(int part) const {
		ConstIterator ptr = data_.cbegin() + part_.end(part), end = this->cend();
		while (ptr < end && ptr[-1].key() == ptr[0].key())
			++ptr;
		return ptr;
	}

	int64_t parts() const {
		return part_.parts;
	}

	const T& front() const {
		return data_.front();
	}

	const T& back() const {
		return data_.back();
	}

	void sort() {
		data_.sort();
	}

private:

	const size_t size_;
	Container data_;
	const Partition<int64_t> part_;

};
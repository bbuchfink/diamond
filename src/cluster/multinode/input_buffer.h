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
#include <type_traits>
#include "util/memory/memory_resource.h"
#include "volume.h"
#include "util/algo/partition.h"
#include "util/parallel/simple_thread_pool.h"

template<typename T>
struct InputBuffer {

	using Container = typename std::conditional<T::POD, std::vector<T>, std::pmr::list<T>>::type;
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
				File in(f[v].path, "rb", File::Flags::DETECT_COMPRESSION);
				in.read((void*) &data_[f[v].oid_begin], f[v].record_count * sizeof(T));
				in.close();
			}
			};
		for (int i = 0; i < std::min(config.threads_, (int)f.size()); ++i)
			pool.spawn(worker);
		pool.join_all();
	}

	InputBuffer(const VolumedFile& f, std::pmr::memory_resource& mem_pool, int parts = config.threads_) :
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
				File in(f[v].path, "rb", File::Flags::DETECT_COMPRESSION);
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
		sort(data_);
	}

	void sort(std::pmr::list<T>& lst) {
		lst.sort();
	}

private:

	const size_t size_;
	Container data_;
	const Partition<int64_t> part_;

};
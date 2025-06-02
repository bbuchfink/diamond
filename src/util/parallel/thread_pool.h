/****
DIAMOND protein aligner
Copyright (C) 2021-2024 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink

Code developed by Benjamin Buchfink <buchfink@gmail.com>

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

#pragma once
#include <array>
#include <thread>
#include <atomic>
#include <vector>
#include <queue>
#include <condition_variable>
#include <functional>
#include "../log_stream.h"

namespace Util { namespace Parallel {

template<typename F, typename... Args>
void pool_worker(std::atomic<size_t> *partition, size_t thread_id, size_t partition_count, F f, Args... args) {
	size_t p;
	while ((p = partition->fetch_add(1, std::memory_order_relaxed)) < partition_count) {
		f(p, thread_id, args...);
	}
}

template<typename F, typename... Args>
void scheduled_thread_pool(size_t thread_count, F f, Args... args) {
	std::atomic<size_t> partition(0);
	std::vector<std::thread> threads;
	for (size_t i = 0; i < thread_count; ++i)
		threads.emplace_back(f, &partition, i, args...);
	for (std::thread &t : threads)
		t.join();
}

template<typename F, typename... Args>
void scheduled_thread_pool_auto(size_t thread_count, size_t partition_count, F f, Args... args) {
	scheduled_thread_pool(thread_count, pool_worker<F, Args...>, partition_count, f, args...);
}

template<typename F>
void launch_threads(int thread_count, F& f) {
	std::vector<std::thread> threads;
	for (int i = 0; i < thread_count; ++i)
		threads.emplace_back(f);
	for (std::thread& t : threads)
		t.join();
}

}}

struct ThreadPool {

	enum { PRIORITY_COUNT = 2 };

	struct TaskSet {
		TaskSet(ThreadPool& thread_pool, int priority):
			priority(priority),
			total_(0),
			finished_(0),
			thread_pool(&thread_pool)
		{}
		void finish() {
			std::unique_lock<std::mutex> lock(thread_pool->mtx_);
			++finished_;
			if (finished())
				cv_.notify_all();
		}
		bool finished() const {
			return total_ == finished_;
		}
		int64_t total() const {
			return total_;
		}
		void run() {
			{
				std::unique_lock<std::mutex> lock(thread_pool->mtx_);
				if (finished())
					return;
			}
			thread_pool->run_set(this);
		}
		template<class F, class... Args>
		void enqueue(F&& f, Args&&... args) {
			thread_pool->enqueue(*this, std::forward<F>(f), std::forward<Args>(args)...);
		}
		const int priority;
	private:		
		int64_t total_, finished_;
		ThreadPool* thread_pool;
		std::condition_variable cv_;
		friend struct ThreadPool;
	};

	struct Task {
		Task():
			task_set(nullptr)
		{}
		Task(std::function<void()> f):
			f(f),
			task_set(nullptr)
		{}
		Task(std::function<void()> f, TaskSet& task_set):
			f(f),
			task_set(&task_set)
		{}
		operator bool() const {
			return f.operator bool();
		}
		std::function<void()> f;
		TaskSet* task_set;
	};

	template<class F, class... Args>
	void enqueue(TaskSet& task_set, F&& f, Args&&... args)
	{
		if (pop_before_enqueue_) {
			for (;;) {
				Task t;
				if (!queue_empty()) {
					std::unique_lock<std::mutex> lock(this->mtx_);
					t = pop_task();
				}
				if (t) {
					t.f();
					if (t.task_set)
						t.task_set->finish();
				}
				else
					break;
			}
		}
		auto task = std::bind(std::forward<F>(f), std::forward<Args>(args)...);
		{
			std::unique_lock<std::mutex> lock(mtx_);
			++task_set.total_;
			tasks_[task_set.priority].emplace([task]() { task(); }, task_set);
			task_set.cv_.notify_one();
		}
	}

	void run_set(TaskSet* task_set) {
		for (;;)
		{
			Task task;

			if (!task_set) {
				if (default_finished_.load(std::memory_order_relaxed) >= default_count_) {
					++threads_finished_;
					return;
				}
				if (!queue_empty()) {
					std::unique_lock<std::mutex> lock(this->mtx_);
					task = pop_task();
				}
				if (!task) {
					const int64_t next = default_begin_.fetch_add(1, std::memory_order_relaxed);
					if (next < default_end_) {
						default_task_(*this, next);
						default_finished_.fetch_add(1, std::memory_order_relaxed);
					}
					continue;
				}
			}
			else {
				std::unique_lock<std::mutex> lock(this->mtx_);
				task_set->cv_.wait(lock,
					[this, task_set] { return !queue_empty(task_set->priority) || task_set->finished(); });
				if (task_set->finished())
					return;
				task = pop_task(task_set->priority);
			}

			task.f();
			if (task.task_set)
				task.task_set->finish();
		}
	}

	ThreadPool(const std::function<void(ThreadPool&, int64_t)>& default_task = std::function<void(ThreadPool&, int64_t)>(), int64_t default_begin = 0, int64_t default_end = 0, bool pop_before_enqueue = false) :
		pop_before_enqueue_(pop_before_enqueue),
		default_end_(default_end),
		default_count_(default_end - default_begin),
		default_task_(default_task),
		mtx_(),
		default_begin_(default_begin),
		default_finished_(0),
		threads_finished_(0)
	{
	}

	void run(int threads, bool heartbeat = false, TaskSet* task_set = nullptr) {
		for (int i = 0; i < threads; ++i)
			workers_.emplace_back([this, task_set] { this->run_set(task_set); });
		if (heartbeat)
			heartbeat_ = std::thread([&]() {
			while (default_finished_ < default_count_) {
				log_stream << "Workers=" << workers_.size() << '/' << threads_finished_ << " begin = " << default_begin_ << " finished = "
					<< default_finished_ << " queue=" << queue_len(0) << '/' << queue_len(1) << std::endl;
				std::this_thread::sleep_for(std::chrono::seconds(1));
			}});
	}

	~ThreadPool()
	{
		join();
	}

	void join() {
		for (std::thread &worker : workers_)
			worker.join();
		workers_.clear();
		if(heartbeat_.joinable())
			heartbeat_.join();
	}

	int64_t queue_len(int priority) const {
		return tasks_[priority].size();
	}

private:

	bool queue_empty(int priority = PRIORITY_COUNT - 1) {
		for (int i = 0; i <= priority; ++i)
			if (!tasks_[i].empty())
				return false;
		return true;
	}

	Task pop_task(int priority = PRIORITY_COUNT - 1) {
		Task task;
		for (int i = 0; i <= priority; ++i)
			if (!tasks_[i].empty()) {
				task = std::move(tasks_[i].front());
				tasks_[i].pop();
				break;
			}
		return task;
	}

	const bool pop_before_enqueue_;
	const int64_t default_end_, default_count_;
	std::array<std::queue<Task>, PRIORITY_COUNT> tasks_;
	std::function<void(ThreadPool&, int64_t)> default_task_;
	std::vector<std::thread> workers_;
	std::thread heartbeat_;
	std::mutex mtx_;
	std::atomic<int64_t> default_begin_, default_finished_, threads_finished_;

};

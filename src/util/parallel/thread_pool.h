/****
DIAMOND protein aligner
Copyright (C) 2021-2022 Max Planck Society for the Advancement of Science e.V.

Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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
#include <numeric>
#include <functional>
#include "../log_stream.h"

namespace Util { namespace Parallel {

template<typename F, typename... Args>
void pool_worker(std::atomic<size_t> *partition, size_t thread_id, size_t partition_count, F f, Args... args) {
	size_t p;
	while ((p = (*partition)++) < partition_count) {
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
			++finished_;
			if (finished()) {
				cv_.notify_all();
				thread_pool->cv_.notify_all();
			}
		}
		bool finished() const {
			return total_ == finished_;
		}
		int64_t total() const {
			return total_;
		}
		void add() {
			++total_;
		}
		void wait() {
			std::unique_lock<std::mutex> lock(mtx_);
			cv_.wait(lock, [this] {return this->finished(); });
		}
		void run() {
			if (finished())
				return;
			thread_pool->run_set(this);
		}
		template<class F, class... Args>
		void enqueue(F&& f, Args&&... args) {
			thread_pool->enqueue(*this, std::forward<F>(f), std::forward<Args>(args)...);
		}
		const int priority;
	private:		
		std::atomic<int64_t> total_, finished_;
		ThreadPool* thread_pool;
		std::mutex mtx_;
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
		auto task = std::bind(std::forward<F>(f), std::forward<Args>(args)...);
		{
			std::unique_lock<std::mutex> lock(mtx_);

			/*if (stop_)
				throw std::runtime_error("enqueue on stopped ThreadPool");*/

			task_set.add();
			tasks_[task_set.priority].emplace([task]() { task(); }, task_set);
		}
		cv_.notify_one();
	}

	void run_set(TaskSet* task_set) {
		for (;;)
		{
			Task task;

			if (!task_set && run_default_) {
				if (!queue_empty()) {
					std::unique_lock<std::mutex> lock(this->mtx_);
					task = pop_task();
				}
				if (!task) {
					++default_started_;
					if (!default_task_(*this))
						run_default_ = false;
					++default_finished_;
					if (!run_default_ && default_started_ == default_finished_) {
						stop_ = true;
						cv_.notify_all();
					}
					continue;
				}
			}
			else {
				std::unique_lock<std::mutex> lock(this->mtx_);
				this->cv_.wait(lock,
					[this, task_set] { return (stop_ && !task_set) || !queue_empty(task_set ? task_set->priority : PRIORITY_COUNT-1) || (task_set && task_set->finished()); });
				if ((stop_ && queue_empty() && !task_set) || (task_set && task_set->finished())) {
					if (!task_set)
						++threads_finished_;
					return;
				}
				task = pop_task(task_set ? task_set->priority : PRIORITY_COUNT - 1);
			}

			task.f();
			if (task.task_set)
				task.task_set->finish();
		}
	}

	ThreadPool(const std::function<bool(ThreadPool&)>& default_task = std::function<bool(ThreadPool&)>()) :
		default_task_(default_task),
		mtx_(),
		stop_(false),
		run_default_(default_task.operator bool()),
		default_started_(0),
		default_finished_(0),
		threads_finished_(0)
	{
	}

	void run(size_t threads, bool heartbeat = false) {
		for (size_t i = 0; i < threads; ++i)
			workers_.emplace_back([this] {	this->run_set(nullptr); });
		if (heartbeat)
			heartbeat_ = std::thread([&]() {
			while (!stop_) {
				log_stream << "Workers=" << workers_.size() << '/' << threads_finished_ << " started = " << default_started_ << " finished = "
					<< default_finished_ << " queue=" << queue_len(0) << '/' << queue_len(1) << std::endl;
				std::this_thread::sleep_for(std::chrono::seconds(1));
			}});
	}

	~ThreadPool()
	{
		{
			std::unique_lock<std::mutex> lock(mtx_);
			stop_ = true;
		}
		cv_.notify_all();
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

	std::array<std::queue<Task>, PRIORITY_COUNT> tasks_;
	std::function<bool(ThreadPool&)> default_task_;
	std::vector<std::thread> workers_;
	std::thread heartbeat_;
	std::mutex mtx_;
	std::condition_variable cv_;
	bool stop_, run_default_;
	std::atomic<int64_t> default_started_, default_finished_, threads_finished_;

};

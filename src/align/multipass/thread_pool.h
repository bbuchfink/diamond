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
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <vector>
#include <functional>
#include <stdexcept>
#include <exception>
#include <iostream>
#include <cstdlib>
#include <utility>

namespace MultipassPipeline {

// Thread pool serving a task queue of generic callables (the query and target tasks, each
// of which chains its successor). The workers wait for tasks until the pool is stopped, which
// is up to the owner as the pool cannot tell whether further tasks will be submitted. Once
// stopped, the workers complete the tasks still queued and exit.
class ThreadPool {
public:

	explicit ThreadPool(unsigned n) {
		if (n == 0) n = 1;
		workers_.reserve(n);
		for (unsigned i = 0; i < n; ++i)
			workers_.emplace_back(&ThreadPool::worker_loop, this, (int)i);
	}

	// Index of the calling worker in [0, n), which a task uses to access data kept per worker without a lock.
	// -1 if the calling thread is not a worker.
	static int worker_index() {
		return worker_index_();
	}

	void stop() {
		{ std::lock_guard<std::mutex> lk(mtx_); stopping_ = true; }
		cv_.notify_all();
	}

	// Waits for the workers to exit, which requires stop() to be called. Not to be called by a worker.
	void join() {
		for (auto& t : workers_)
			if (t.joinable()) t.join();
	}

	~ThreadPool() {
		join();
	}

	ThreadPool(const ThreadPool&) = delete;
	ThreadPool& operator=(const ThreadPool&) = delete;

	void submit(std::function<void()> task) {         // callable from any thread, incl. workers
		{
			std::lock_guard<std::mutex> lk(mtx_);
			if (stopping_) throw std::runtime_error("submit on stopped pool");
			task_queue_.push(std::move(task));
		}
		cv_.notify_one();
	}

private:

	static int& worker_index_() {
		static thread_local int index = -1;
		return index;
	}

	void worker_loop(int index) {
		worker_index_() = index;
		for (;;) {
			std::function<void()> task;
			{
				std::unique_lock<std::mutex> lk(mtx_);
				cv_.wait(lk, [this] { return stopping_ || !task_queue_.empty(); });
				if (task_queue_.empty())
					return;                            // stopping and nothing left to do
				task = std::move(task_queue_.front());
				task_queue_.pop();
			}
			try {
				task();                                // run outside the lock
			}
			catch (std::exception& e) {
				std::cerr << "Error: " << e.what() << std::endl;
				exit(EXIT_FAILURE);
			}
		}
	}

	std::mutex mtx_;                                   // lock of the task queue, guards task_queue_ and stopping_
	std::condition_variable cv_;
	std::queue<std::function<void()>> task_queue_;
	std::vector<std::thread> workers_;
	bool stopping_ = false;
};

}

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
#include <utility>

namespace ExtensionPipeline {

// Thread pool serving two work queues: a task queue of generic callables (the query
// and target tasks, each of which chains its successor) and the extension queue
// holding the individual extensions submitted by the target tasks. A worker takes an
// extension whenever one is available and falls back to the task queue otherwise, so
// new extensions are only generated once the pending ones have been handed out.
template<typename Ext>
class ThreadPool {
public:

    ThreadPool(unsigned n, std::function<void(Ext&)> extension_fn) :
        extension_fn_(std::move(extension_fn))
    {
        if (n == 0) n = 1;
        workers_.reserve(n);
        for (unsigned i = 0; i < n; ++i)
            workers_.emplace_back(&ThreadPool::worker_loop, this);
    }

    void stop() {
        { std::lock_guard<std::mutex> lk(mu_); stopping_ = true; }
        cv_.notify_all();
    }

    ~ThreadPool() {
        for (auto& t : workers_) t.join();
    }

    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

    void submit(std::function<void()> task) {         // callable from any thread, incl. workers
        {
            std::lock_guard<std::mutex> lk(mu_);
            if (stopping_) throw std::runtime_error("submit on stopped pool");
            task_queue_.push(std::move(task));
        }
        cv_.notify_one();
    }

    void submit_extension(Ext extension) {            // callable from any thread, incl. workers
        {
            std::lock_guard<std::mutex> lk(mu_);
            if (stopping_) throw std::runtime_error("submit on stopped pool");
            extension_queue_.push(std::move(extension));
        }
        cv_.notify_one();
    }

private:

    void worker_loop() {
        for (;;) {
            Ext extension;
            std::function<void()> task;
            bool have_extension = false;
            {
                std::unique_lock<std::mutex> lk(mu_);
                cv_.wait(lk, [this] { return stopping_ || !extension_queue_.empty() || !task_queue_.empty(); });
                if (!extension_queue_.empty()) {
                    extension = std::move(extension_queue_.front());
                    extension_queue_.pop();
                    have_extension = true;
                }
                else if (!task_queue_.empty()) {
                    task = std::move(task_queue_.front());
                    task_queue_.pop();
                }
                else
                    return;                            // stopping and nothing left to do
            }
            if (have_extension)                        // run outside the lock
                extension_fn_(extension);
            else
                task();
        }
    }

    std::function<void(Ext&)> extension_fn_;
    std::mutex mu_;
    std::condition_variable cv_;
    std::queue<std::function<void()>> task_queue_;
    std::queue<Ext> extension_queue_;
    std::vector<std::thread> workers_;
    bool stopping_ = false;
};

}

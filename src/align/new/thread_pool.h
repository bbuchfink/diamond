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
#include <cstdint>
#include "extension.h"

namespace ExtensionPipeline {

// Thread pool serving two work queues: a task queue of generic callables (the query
// and target tasks, each of which chains its successor) and the extension queue
// holding the individual extensions submitted by the target tasks. Extensions are
// computed in batches: a worker calls the batch function, which pops all queued
// extensions at once, as soon as at least MIN_BATCH of them are available, and runs a
// task otherwise. Smaller batches are only computed once no task is queued or running
// any more, i.e. no further extensions can be generated.
// The mutex of the extension queue guards all state of the pool, so the batch function
// must not be called with it held.
// Two workers may decide to compute a batch before the first one has popped the queue,
// so the batch function has to cope with an empty queue.
class ThreadPool {
public:

    enum { MIN_BATCH = 16 };

    ThreadPool(unsigned n, std::function<void(ExtensionQueue&)> batch_fn) :
        batch_fn_(std::move(batch_fn))
    {
        if (n == 0) n = 1;
        workers_.reserve(n);
        for (unsigned i = 0; i < n; ++i)
            workers_.emplace_back(&ThreadPool::worker_loop, this);
    }

    void stop() {
        { std::lock_guard<std::mutex> lk(extension_queue_.mtx); stopping_ = true; }
        cv_.notify_all();
    }

    ~ThreadPool() {
        for (auto& t : workers_) t.join();
    }

    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

    void submit(std::function<void()> task) {         // callable from any thread, incl. workers
        {
            std::lock_guard<std::mutex> lk(extension_queue_.mtx);
            if (stopping_) throw std::runtime_error("submit on stopped pool");
            task_queue_.push(std::move(task));
        }
        cv_.notify_one();
    }

    void submit_extension(Extension&& extension) {    // callable from any thread, incl. workers
        bool notify;
        {
            std::lock_guard<std::mutex> lk(extension_queue_.mtx);
            if (stopping_) throw std::runtime_error("submit on stopped pool");
            extension_queue_.extensions.push_back(std::move(extension));
            notify = extension_queue_.extensions.size() >= MIN_BATCH || (task_queue_.empty() && running_tasks_ == 0);
        }
        if (notify)
            cv_.notify_one();
    }

private:

    // Requires the lock to be held.
    bool batch_ready() const {
        const std::vector<Extension>& extensions = extension_queue_.extensions;
        return !extensions.empty()
            && (extensions.size() >= MIN_BATCH || (task_queue_.empty() && running_tasks_ == 0) || stopping_);
    }

    void worker_loop() {
        for (;;) {
            std::function<void()> task;
            bool have_batch = false;
            {
                std::unique_lock<std::mutex> lk(extension_queue_.mtx);
                cv_.wait(lk, [this] { return stopping_ || batch_ready() || !task_queue_.empty(); });
                if (batch_ready())
                    have_batch = true;
                else if (!task_queue_.empty()) {
                    task = std::move(task_queue_.front());
                    task_queue_.pop();
                    ++running_tasks_;
                }
                else
                    return;                            // stopping and nothing left to do
            }
            if (have_batch)                            // pops the queue itself, run outside the lock
                batch_fn_(extension_queue_);
            else {
                task();
                bool notify;
                {
                    std::lock_guard<std::mutex> lk(extension_queue_.mtx);
                    --running_tasks_;
                    notify = batch_ready();
                }
                if (notify)
                    cv_.notify_one();
            }
        }
    }

    std::function<void(ExtensionQueue&)> batch_fn_;
    ExtensionQueue extension_queue_;
    std::condition_variable cv_;
    std::queue<std::function<void()>> task_queue_;
    std::vector<std::thread> workers_;
    int64_t running_tasks_ = 0;                        // tasks taken from the queue and not yet finished
    bool stopping_ = false;
};

}
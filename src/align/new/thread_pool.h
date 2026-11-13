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
#include <algorithm>
#include <iterator>
#include <condition_variable>
#include <queue>
#include <vector>
#include <deque>
#include <array>
#include <functional>
#include <stdexcept>
#include <utility>
#include <cstdint>
#include <atomic>
#include "extension.h"

namespace ExtensionPipeline {

// Thread pool serving a task queue of generic callables (the query and target tasks, each
// of which chains its successor) and the extension queues holding the extensions submitted
// by the target tasks one at a time, partitioned by band (see ExtensionQueue). Extensions are
// computed in batches by batch tasks, which run through the task queue like any other task.
// Whenever a task has submitted extensions, the worker running it queues batch tasks for the extension queues it
// submitted to, until the queued batch tasks of a queue account for min_batch extensions each. A batch task
// pops min_batch extensions from its queue and calls the batch function with them and the queue.
// The batch function may keep popping extensions submitted to its queue while it runs, so a batch task may
// find its queue holding less than min_batch extensions, in which case it does nothing.
// Smaller batches are only computed once no task is queued or running any more, i.e. no further extensions
// can be generated: the worker finishing the last task then queues a flush task for every non-empty queue, which
// pops the remaining extensions regardless of their number.
// Submitting and popping an extension only takes the lock of one shard of its extension queue, and the batch tasks of
// a queue are counted without a lock. Only queueing a task takes the lock of the task queue. The batch tasks are
// scheduled from the size of the queue, which changes concurrently, so a queue may hold min_batch extensions or more
// without a batch task accounting for them. This only delays them: they are computed by later batch tasks, by
// anchored swipes popping from the queue, or at the latest by the flush tasks.
class ThreadPool {
public:

    using BatchFn = std::function<void(std::vector<Extension>&, ExtensionQueue&)>;

    // band_max has to bound the band of any extension submitted (see queue_band).
    ThreadPool(unsigned n, int64_t min_batch, Loc band_max, BatchFn batch_fn) :
        batch_fn_(std::move(batch_fn)),
        min_batch_((size_t)std::max(min_batch, (int64_t)1))
    {
        if (n == 0) n = 1;
        // One shard per worker, so that a worker submitting extensions has a shard of its own.
        for (int k = 0; k < QUEUE_COUNT; ++k) {
            extension_queues_.emplace_back(queue_band(k, band_max), (int)n);
            scheduled_[k].store(0, std::memory_order_relaxed);
        }
        workers_.reserve(n);
        for (unsigned i = 0; i < n; ++i)
            workers_.emplace_back(&ThreadPool::worker_loop, this, (int)i);
    }

    void stop() {
        { std::lock_guard<std::mutex> lk(mtx_); stopping_ = true; }
        cv_.notify_all();
    }

    ~ThreadPool() {
        for (auto& t : workers_) t.join();
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

    // Pushes the extension to the extension queue of its band, which makes it available to the anchored swipes
    // computing the queue right away. The batch tasks of the queue are topped up when the task of the calling worker
    // is done (see worker_loop), or right here if the caller is not a worker.
    void submit_extension(Extension&& extension) {    // callable from any thread, incl. workers
        if (stopping_.load(std::memory_order_relaxed)) throw std::runtime_error("submit on stopped pool");
        const int k = queue_index(extension.diag_spread());
        extension_queues_[k].push(std::move(extension));
        if (Submitted* const submitted = worker_submitted())
            (*submitted)[k] = true;
        else
            schedule_batches(k);
    }

private:

    // Extension queues the task run by a worker has submitted to.
    using Submitted = std::array<bool, QUEUE_COUNT>;

    // The array above of the calling thread if it is a worker, nullptr otherwise.
    static Submitted*& worker_submitted() {
        static thread_local Submitted* submitted = nullptr;
        return submitted;
    }

    // Wakes up workers for the given number of tasks just queued.
    void notify(size_t tasks) {
        if (tasks == 1)
            cv_.notify_one();
        else if (tasks > 1)
            cv_.notify_all();
    }

    // Queues batch tasks for extension queue k until they account for the extensions it holds, i.e. the queue
    // holds less than min_batch extensions more than them.
    void schedule_batches(int k) {
        const size_t n = extension_queues_[k].size();
        size_t scheduled = scheduled_[k].load(std::memory_order_relaxed), tasks = 0;
        // Counted as scheduled before the tasks are queued, which is safe as a batch task cannot start before.
        while (n >= (scheduled + 1) * min_batch_)
            if (scheduled_[k].compare_exchange_weak(scheduled, scheduled + 1, std::memory_order_relaxed)) {
                ++scheduled;
                ++tasks;
            }
        if (tasks == 0)
            return;
        {
            std::lock_guard<std::mutex> lk(mtx_);
            for (size_t i = 0; i < tasks; ++i)
                task_queue_.push([this, k] { run_batch(k, false); });
        }
        notify(tasks);
    }

    // Requires the lock of the task queue to be held and no task to be queued or running. Queues a flush task for
    // every non-empty extension queue and returns their number.
    // The sizes of the queues are exact, as no thread can change them while no task is running, and every change is
    // visible: the task having made it finished by decrementing running_tasks_ under the lock of the task queue, which
    // the caller has acquired since.
    size_t schedule_flushes() {
        size_t tasks = 0;
        for (int k = 0; k < QUEUE_COUNT; ++k)
            if (extension_queues_[k].size() > 0) {
                task_queue_.push([this, k] { run_batch(k, true); });
                ++tasks;
            }
        return tasks;
    }

    // Batch task (flush task if `flush`) of extension queue k. Pops min_batch extensions if the queue holds at least
    // that many, or whatever extensions it holds for a flush task, and calls the batch function with them.
    // A flush task pops at most min_batch extensions as well. As the batch function keeps popping from the queue until
    // it is empty, this empties the queue. If not, the worker finishing the flush task queues another one.
    void run_batch(int k, bool flush) {
        ExtensionQueue& queue = extension_queues_[k];
        if (!flush) {
            scheduled_[k].fetch_sub(1, std::memory_order_relaxed);
            if (queue.size() < min_batch_)
                return;
        }
        std::vector<Extension> batch(min_batch_);
        int shard = ExtensionQueue::thread_shard();
        batch.resize(queue.pop(shard, batch.data(), min_batch_));
        if (batch.empty())
            return;
        batch_fn_(batch, queue);
    }

    void worker_loop(int index) {
        Submitted submitted{};
        worker_submitted() = &submitted;
        ExtensionQueue::thread_shard() = index;
        for (;;) {
            std::function<void()> task;
            {
                std::unique_lock<std::mutex> lk(mtx_);
                cv_.wait(lk, [this] { return stopping_ || !task_queue_.empty(); });
                if (task_queue_.empty())
                    return;                            // stopping and nothing left to do
                task = std::move(task_queue_.front());
                task_queue_.pop();
                ++running_tasks_;
            }
            task();                                    // run outside the lock
            // Before the task counts as finished, so that the batch tasks are queued when the flush tasks are considered.
            for (int k = 0; k < QUEUE_COUNT; ++k)
                if (submitted[k]) {
                    submitted[k] = false;
                    schedule_batches(k);
                }
            size_t tasks = 0;
            {
                std::lock_guard<std::mutex> lk(mtx_);
                if (--running_tasks_ == 0 && task_queue_.empty())
                    tasks = schedule_flushes();
            }
            notify(tasks);
        }
    }

    BatchFn batch_fn_;
    const size_t min_batch_;                           // extensions popped per batch
    std::mutex mtx_;                                   // lock of the task queue, guards task_queue_, running_tasks_ and stopping_
    std::deque<ExtensionQueue> extension_queues_;      // indexed by queue_index (a deque as the queues are not movable)
    // Batch tasks queued but not yet started, per extension queue.
    std::array<std::atomic<size_t>, QUEUE_COUNT> scheduled_;
    std::condition_variable cv_;
    std::queue<std::function<void()>> task_queue_;
    std::vector<std::thread> workers_;
    int64_t running_tasks_ = 0;                        // tasks (incl. batch tasks) taken from the queue and not yet finished
    // Written under the lock of the task queue, so that waiting workers are not missed, and read without it by submit_extension.
    std::atomic<bool> stopping_{ false };
};

}

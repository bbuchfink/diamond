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

#include <thread>
#include <random>
#include <vector>
#include <atomic>
#include <iostream>
#include <cstdint>
#include <sstream>
#include "../util/data_structures/queue.h"
#include "util/parallel/filestack.h"
#include "basic/config.h"

using std::ostringstream;
using std::unique_ptr;
using std::thread;
using std::endl;
using std::vector;

struct QueueStressTestResult {
    bool passed;
    size_t items_sent;
    size_t items_received;
    uint64_t expected_checksum;
    uint64_t received_checksum;
};

static QueueStressTestResult test_many_producers_one_consumer(int thread_count, size_t items_per_producer) {
    const int producer_count = thread_count - 1;
    const int consumer_count = 1;
    const int64_t poison_pill = -1;
    const size_t queue_capacity = 1024;

    if (producer_count < 1) {
        return { false, 0, 0, 0, 0 };
    }

    Queue<int64_t> queue(queue_capacity, producer_count, consumer_count, poison_pill);

    std::atomic<size_t> total_sent(0);
    std::atomic<size_t> total_sent2(0);
    std::atomic<size_t> total_received(0);
    std::atomic<uint64_t> sent_checksum(0);
    std::atomic<uint64_t> received_checksum(0);

    std::vector<std::thread> producers;
    producers.reserve(producer_count);

    for (int p = 0; p < producer_count; ++p) {
        producers.emplace_back([&queue, &total_sent, &total_sent2, &sent_checksum, p, items_per_producer]() {
            uint64_t local_checksum = 0;
            for (size_t i = 0; i < items_per_producer; ++i) {
                int64_t value = static_cast<int64_t>(p) * static_cast<int64_t>(items_per_producer) + static_cast<int64_t>(i);
                queue.enqueue(value);
                local_checksum += static_cast<uint64_t>(value);
                size_t x = total_sent2.fetch_add(1, std::memory_order_relaxed);
            }
            total_sent.fetch_add(items_per_producer, std::memory_order_relaxed);
            sent_checksum.fetch_add(local_checksum, std::memory_order_relaxed);
            queue.enqueue(-1);            
        });
    }

    std::thread consumer([&queue, &total_received, &received_checksum]() {
        int64_t value;
        uint64_t local_checksum = 0;
        size_t count = 0;
        while (queue.wait_and_dequeue(value)) {
            local_checksum += static_cast<uint64_t>(value);
            ++count;
        }
        total_received.fetch_add(count, std::memory_order_relaxed);
        received_checksum.fetch_add(local_checksum, std::memory_order_relaxed);
    });

    int j = 0;
    for (auto& t : producers) {
        t.join();
    }
    consumer.join();

    const size_t expected_count = static_cast<size_t>(producer_count) * items_per_producer;
    const size_t sent = total_sent.load();
    const size_t received = total_received.load();
    const uint64_t sent_sum = sent_checksum.load();
    const uint64_t received_sum = received_checksum.load();

    bool passed = (sent == expected_count) && (received == expected_count) && (sent_sum == received_sum);

    return { passed, sent, received, sent_sum, received_sum };
}

static QueueStressTestResult test_one_producer_many_consumers(int thread_count, size_t total_items) {
    const int producer_count = 1;
    const int consumer_count = thread_count - 1;
    const int64_t poison_pill = -1;
    const size_t queue_capacity = 1024;

    if (consumer_count < 1) {
        return { false, 0, 0, 0, 0 };
    }

    Queue<int64_t> queue(queue_capacity, producer_count, consumer_count, poison_pill);

    std::atomic<size_t> total_received(0);
    std::atomic<uint64_t> received_checksum(0);
    uint64_t expected_checksum = 0;

    std::thread producer([&queue, total_items, &expected_checksum]() {
        for (size_t i = 0; i < total_items; ++i) {
            int64_t value = static_cast<int64_t>(i);
            queue.enqueue(value);
            expected_checksum += static_cast<uint64_t>(value);
        }
        queue.close();
    });

    std::vector<std::thread> consumers;
    consumers.reserve(consumer_count);

    for (int c = 0; c < consumer_count; ++c) {
        consumers.emplace_back([&queue, &total_received, &received_checksum]() {
            int64_t value;
            uint64_t local_checksum = 0;
            size_t count = 0;
            while (queue.wait_and_dequeue(value)) {
                local_checksum += static_cast<uint64_t>(value);
                ++count;
            }
            total_received.fetch_add(count, std::memory_order_relaxed);
            received_checksum.fetch_add(local_checksum, std::memory_order_relaxed);
        });
    }

    producer.join();
    for (auto& t : consumers) {
        t.join();
    }

    const size_t received = total_received.load();
    const uint64_t received_sum = received_checksum.load();

    bool passed = (received == total_items) && (expected_checksum == received_sum);

    return { passed, total_items, received, expected_checksum, received_sum };
}

int run_queue_stress_test() {
    const int thread_count = static_cast<int>(std::thread::hardware_concurrency());
    const size_t items_per_producer = 300;
    const size_t total_items = 10000;

    std::cout << "Queue Stress Test" << std::endl;
    std::cout << "=================" << std::endl;
    std::cout << "Hardware threads: " << thread_count << std::endl;
    std::cout << std::endl;

    if (thread_count < 2) {
        std::cout << "Error: Need at least 2 threads for stress test" << std::endl;
        return 1;
    }

    int failures = 0;

    // Test 1: Many producers, one consumer
    {
        std::cout << "Test 1: Many producers (" << (thread_count - 1) << "), one consumer" << std::endl;
        std::cout << "  Items per producer: " << items_per_producer << std::endl;
        std::cout << "  Total items: " << (thread_count - 1) * items_per_producer << std::endl;

        auto result = test_many_producers_one_consumer(thread_count, items_per_producer);

        std::cout << "  Items sent: " << result.items_sent << std::endl;
        std::cout << "  Items received: " << result.items_received << std::endl;
        std::cout << "  Expected checksum: " << result.expected_checksum << std::endl;
        std::cout << "  Received checksum: " << result.received_checksum << std::endl;
        std::cout << "  Result: " << (result.passed ? "PASSED" : "FAILED") << std::endl;
        std::cout << std::endl;

        if (!result.passed) ++failures;
    }

    // Test 2: One producer, many consumers
    {
        std::cout << "Test 2: One producer, many consumers (" << (thread_count - 1) << ")" << std::endl;
        std::cout << "  Total items: " << total_items << std::endl;

        auto result = test_one_producer_many_consumers(thread_count, total_items);

        std::cout << "  Items sent: " << result.items_sent << std::endl;
        std::cout << "  Items received: " << result.items_received << std::endl;
        std::cout << "  Expected checksum: " << result.expected_checksum << std::endl;
        std::cout << "  Received checksum: " << result.received_checksum << std::endl;
        std::cout << "  Result: " << (result.passed ? "PASSED" : "FAILED") << std::endl;
        std::cout << std::endl;

        if (!result.passed) ++failures;
    }

    std::cout << "=================" << std::endl;
    std::cout << "Tests passed: " << (2 - failures) << "/2" << std::endl;

    return failures;
}

void filestack() {
    unique_ptr<FileStack> stack(new FileStack("test.tsv"));
    vector<thread> workers;
    auto worker = [&](int thread_id) {
        for (int i = 0; i < 100; ++i) {
            std::random_device r;
            std::default_random_engine e1(r());
            std::uniform_int_distribution<int> uniform_dist(1, 999999999);
            ostringstream ss;
            ss << uniform_dist(e1) << '\t' << uniform_dist(e1) << '\t' << uniform_dist(e1) << '\t' << uniform_dist(e1) << '\n';
            stack->push(ss.str());
        }
        };
    for (int i = 0; i < config.threads_; ++i)
        workers.emplace_back(worker, i);
    for (auto& t : workers) {
        t.join();
    }
}
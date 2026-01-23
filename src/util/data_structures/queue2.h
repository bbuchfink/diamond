#pragma once
#include <mutex>
#include <condition_variable>
#include <deque>
#include <vector>
#include <stdexcept>
#include <type_traits>

template <class T>
class Queue {
public:
    using value_type = T;

    static_assert(std::is_default_constructible_v<T>, "T must be default constructible");

    // Constructor
    explicit Queue(size_t capacity, int producer_count, int consumer_count, T poison_pill)
        : capacity_(capacity ? capacity : 1),
        producer_count_(producer_count),
        consumer_count_(consumer_count),
        poison_pill_(std::move(poison_pill)),
        pills_received_(0) {
    }

    // Default destructor is fine as std::deque handles object destruction,
    // but to match the "drain" behavior of the original if consumers quit early:
    ~Queue() {
        std::lock_guard<std::mutex> lock(mtx_);
        queue_.clear();
    }

    Queue(const Queue&) = delete;
    Queue& operator=(const Queue&) = delete;

    void enqueue(const T& v) {
        std::unique_lock<std::mutex> lock(mtx_);
        // Wait until there is space
        not_full_.wait(lock, [this] { return queue_.size() < capacity_; });

        queue_.push_back(v);

        // Notify a waiting consumer
        not_empty_.notify_one();
    }

    void enqueue(T&& v) {
        std::unique_lock<std::mutex> lock(mtx_);
        not_full_.wait(lock, [this] { return queue_.size() < capacity_; });

        queue_.push_back(std::move(v));

        not_empty_.notify_one();
    }

    bool try_dequeue(T& out) {
        std::lock_guard<std::mutex> lock(mtx_);
        if (queue_.empty()) {
            return false;
        }

        out = std::move(queue_.front());
        queue_.pop_front();

        // Notify a waiting producer
        not_full_.notify_one();
        return true;
    }

    bool wait_and_dequeue(T& out) {
        std::unique_lock<std::mutex> lock(mtx_);

        while (true) {
            // Wait until the queue is not empty
            not_empty_.wait(lock, [this] { return !queue_.empty(); });

            // Retrieve the item (but don't pop immediately until we check logic)
            T val = std::move(queue_.front());
            queue_.pop_front();

            // Notify producers that space is available
            not_full_.notify_one();

            // Check for poison pill
            if (val == poison_pill_) {
                // Scenario: Many Producers, 1 Consumer (1C-MP) logic
                if (producer_count_ > 1) {
                    pills_received_++;

                    if (pills_received_ == static_cast<size_t>(producer_count_)) {
                        return false; // All producers have signed off
                    }

                    // We haven't seen enough pills yet. 
                    // Swallow this pill and continue waiting in the loop.
                    continue;
                }

                // Scenario: 1 Producer or Default
                return false;
            }

            // Real data received
            out = std::move(val);
            return true;
        }
    }

    bool empty() const {
        std::lock_guard<std::mutex> lock(mtx_);
        return queue_.empty();
    }

    size_t capacity() const {
        return capacity_;
    }

    size_t approx_size() const {
        std::lock_guard<std::mutex> lock(mtx_);
        return queue_.size();
    }

    void close() {
        // Enqueue one poison pill for every consumer.
        // Since enqueue is thread-safe and blocking, this ensures
        // we don't overrun the buffer even during shutdown.
        for (int i = 0; i < consumer_count_; ++i) {
            enqueue(poison_pill_);
        }
    }

    int producer_count() const noexcept {
        return producer_count_;
	}

private:
    // Synchronization primitives
    mutable std::mutex mtx_;
    std::condition_variable not_full_;
    std::condition_variable not_empty_;

    // Underlying container
    std::deque<T> queue_;
    const size_t capacity_;

    // Logic configuration
    const int producer_count_;
    const int consumer_count_;
    const T poison_pill_;

    // Internal state (protected by mtx_)
    size_t pills_received_;
};
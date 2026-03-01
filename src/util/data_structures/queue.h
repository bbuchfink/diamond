/****
Copyright (C) 2012-2026 Benjamin J. Buchfink <buchfink@gmail.com>

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
#include <atomic>
#include <assert.h>
#include <type_traits>
#include <new>
#include "../parallel/semaphore.h"

template <class T>
class Queue {

public:
    using value_type = T;

#if _MSC_VER || __GNUC__ >= 5
    static_assert(std::is_default_constructible_v<T>, "T must be default constructible");
#endif

    explicit Queue(size_t capacity, int producer_count, int consumer_count, T poison_pill)
        : capacity_(round_up_to_power_of_two(capacity ? capacity : 1)),
        mask_(capacity_ - 1),
        producer_count_(producer_count),
        consumer_count_(consumer_count),
        poison_pill_(std::move(poison_pill)),
#if defined(__cpp_aligned_new) && __cpp_aligned_new >= 201606L
        buffer_(static_cast<cell_t*>(::operator new[](capacity_ * sizeof(cell_t), std::align_val_t{ alignof(cell_t) }))),
#else
        buffer_(static_cast<cell_t*>(::operator new[](capacity_ * sizeof(cell_t)))),        
#endif
        enqueue_pos_(0),
        dequeue_pos_(0),
        pills_received_(0) {
        assert(producer_count == 1 || consumer_count == 1);
        for (size_t i = 0; i < capacity_; ++i) {
            new (&buffer_[i]) cell_t();
            buffer_[i].seq.store(i, std::memory_order_relaxed);
        }
    }

    ~Queue() {
        T tmp;
        while (try_dequeue(tmp)) {}

        for (size_t i = 0; i < capacity_; ++i) {
            buffer_[i].~cell_t();
        }
#if defined(__cpp_aligned_new) && __cpp_aligned_new >= 201606L
        ::operator delete[](buffer_, std::align_val_t{ alignof(cell_t) });
#else
        ::operator delete[](buffer_);
#endif
    }

    Queue(const Queue&) = delete;
    Queue& operator=(const Queue&) = delete;

    void enqueue(const T& v) { emplace_impl(v); }
    void enqueue(T&& v) { emplace_impl(std::move(v)); }
   
    bool wait_and_dequeue(T& out) {
        while (true) {
            items_.acquire();
            while (!try_dequeue(out)) {
                std::this_thread::yield();
            }
            spaces_.release();
            if (out == poison_pill_) {
                if (producer_count_ > 1) {
                    size_t seen = pills_received_.fetch_add(1, std::memory_order_relaxed) + 1;
                    if (seen == static_cast<size_t>(producer_count_)) {
                        return false;
                    }
                    continue;
                }
                return false;
            }

            return true;
        }
    }

    bool empty() const {
        return approx_size() == 0;
    }

    size_t capacity() const { return capacity_; }

    size_t approx_size() const {
        size_t enq = enqueue_pos_.load(std::memory_order_relaxed);
        size_t deq = dequeue_pos_.load(std::memory_order_relaxed);
        return enq >= deq ? (enq - deq) : 0;
    }

    // call once per producer
    void close() {
        for (int i = 0; i < consumer_count_; ++i) {
            enqueue(poison_pill_);
        }
    }

    int producer_count() const {
		return producer_count_;
    }

private:

    struct cell_t {
        struct Storage {
            alignas(T) unsigned char storage[sizeof(T)];
        };
        std::atomic<size_t> seq;
        Storage storage;
        char pad[(64 - ((sizeof(std::atomic<size_t>) + sizeof(Storage)) % 64)) % 64];
        cell_t() : seq(0) {}
    };

    static constexpr bool is_power_of_two(size_t x) {
        return (x & (x - 1)) == 0;
    }

    static size_t round_up_to_power_of_two(size_t x) {
        if (is_power_of_two(x)) return x;
        --x;
        for (size_t i = 1; i < sizeof(size_t) * 8; i <<= 1)
            x |= x >> i;
        return x + 1;
    }

    static T* cell_value_ptr(cell_t* c) noexcept {
#if defined(__cpp_lib_launder) && __cpp_lib_launder >= 201606L
        return std::launder(reinterpret_cast<T*>(c->storage.storage));
#else
#warning "old compiler missing support for std::launder"
        return reinterpret_cast<T*>(&c->storage);
#endif
    }

    size_t index(size_t pos) const noexcept { return pos & mask_; }

    template <class... Args>
    void emplace_impl(Args&&... args) {        
        spaces_.acquire();
        size_t pos = enqueue_pos_.load(std::memory_order_relaxed);
        for (;;) {
            cell_t* cell = &buffer_[index(pos)];
            size_t seq = cell->seq.load(std::memory_order_acquire);
            intptr_t dif = static_cast<intptr_t>(seq) - static_cast<intptr_t>(pos);

            if (dif == 0) {
                if (enqueue_pos_.compare_exchange_weak(pos, pos + 1, std::memory_order_relaxed)) {
                    new (cell_value_ptr(cell)) T(std::forward<Args>(args)...);
                    cell->seq.store(pos + 1, std::memory_order_release);
                    items_.release();                    
                    return;
                }
            }
            else {
                std::this_thread::yield();
                pos = enqueue_pos_.load(std::memory_order_relaxed);
            }
        }
    }

    bool try_dequeue(T& out) {
        size_t pos = dequeue_pos_.load(std::memory_order_relaxed);
        for (;;) {
            cell_t* cell = &buffer_[index(pos)];
            size_t seq = cell->seq.load(std::memory_order_acquire);
            intptr_t dif = static_cast<intptr_t>(seq) - static_cast<intptr_t>(pos + 1);
            if (dif == 0) {
                if (dequeue_pos_.compare_exchange_weak(
                    pos, pos + 1,
                    std::memory_order_relaxed,
                    std::memory_order_relaxed)) {
                    T* p = cell_value_ptr(cell);
                    out = std::move(*p);
                    p->~T();
                    cell->seq.store(pos + capacity_, std::memory_order_release);
                    return true;
                }
            }
            else if (dif < 0) {
                return false;
            }
            else {
                std::this_thread::yield();
                pos = dequeue_pos_.load(std::memory_order_relaxed);
            }
        }
    }

    const size_t capacity_;
    const size_t mask_;
    const int producer_count_;
    const int consumer_count_;
    const T poison_pill_;
    cell_t* const buffer_;
    alignas(64) std::atomic<size_t> enqueue_pos_;
    char _pad0[(64 - sizeof(enqueue_pos_) % 64) % 64]{};
    alignas(64) std::atomic<size_t> dequeue_pos_;
    char _pad1[(64 - sizeof(dequeue_pos_) % 64) % 64]{};
    CountingSemaphore<> items_{ (ptrdiff_t)0 };
    CountingSemaphore<> spaces_{ (ptrdiff_t)capacity_ };
    std::atomic<size_t> pills_received_;

};
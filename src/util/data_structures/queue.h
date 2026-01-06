/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <atomic>
#include <type_traits>
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
        buffer_(static_cast<cell_t*>(::operator new[](capacity_ * sizeof(cell_t)))),
        enqueue_pos_(0),
        dequeue_pos_(0),
        pills_received_(0) {

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
        ::operator delete[](buffer_);
    }

    Queue(const Queue&) = delete;
    Queue& operator=(const Queue&) = delete;

    void enqueue(const T& v) { emplace_impl(v); }
    void enqueue(T&& v) { emplace_impl(std::move(v)); }

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
                pos = dequeue_pos_.load(std::memory_order_relaxed);
            }
        }
    }

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

    void close() {
        for (int i = 0; i < consumer_count_; ++i) {
            enqueue(poison_pill_);
        }
    }

private:

    struct cell_t {
        struct Storage {
            alignas(T) unsigned char storage[sizeof(T)];
        };
        std::atomic<size_t> seq;
        Storage storage;
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
        return reinterpret_cast<T*>(&c->storage);
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
            else if (dif < 0) {
                pos = enqueue_pos_.load(std::memory_order_relaxed);
            }
            else {
                pos = enqueue_pos_.load(std::memory_order_relaxed);
            }
        }
    }

    const size_t capacity_;
    const size_t mask_;
    const int producer_count_;
    const int consumer_count_;
    const T poison_pill_;
    std::atomic<bool> closed_{ false };
    cell_t* const buffer_;
    alignas(64) std::atomic<size_t> enqueue_pos_;
    char _pad0[64 - sizeof(enqueue_pos_) % 64]{};
    alignas(64) std::atomic<size_t> dequeue_pos_;
    char _pad1[64 - sizeof(dequeue_pos_) % 64]{};
    std::counting_semaphore<> items_{ (ptrdiff_t)0 };
    std::counting_semaphore<> spaces_{ (ptrdiff_t)capacity_ };
    std::atomic<size_t> pills_received_;

};
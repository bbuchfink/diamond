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

#pragma once
#include <atomic>
#include <new>
#include <type_traits>
#include <utility>
#include <semaphore>

template <class T>
class SingleProducerQueue {
public:
    using value_type = T;

    explicit SingleProducerQueue(std::size_t capacity, int consumer_count)
        : capacity_(round_up_to_power_of_two(capacity ? capacity : 1)),
        mask_(capacity_ - 1),
		consumer_count_(consumer_count),        
        buffer_(static_cast<cell_t*>(::operator new[](capacity_ * sizeof(cell_t)))),
        enqueue_pos_(0),
        dequeue_pos_(0) {
        for (std::size_t i = 0; i < capacity_; ++i) {
            new (&buffer_[i]) cell_t();
            buffer_[i].seq.store(i, std::memory_order_relaxed);
        }
    }

    ~SingleProducerQueue() {
        T tmp;
        while (try_dequeue(tmp)) { }

        for (std::size_t i = 0; i < capacity_; ++i) {
            buffer_[i].~cell_t();
        }
        ::operator delete[](buffer_);
    }

    SingleProducerQueue(const SingleProducerQueue&) = delete;
    SingleProducerQueue& operator=(const SingleProducerQueue&) = delete;

    bool try_enqueue(const T& v) { return emplace_impl(v); }
    bool try_enqueue(T&& v) { return emplace_impl(std::move(v)); }

    template <class... Args>
    bool try_emplace(Args&&... args) {
        return emplace_impl(std::forward<Args>(args)...);
    }

    template <class... Args>
    bool emplace(Args&&... args) { return try_emplace(std::forward<Args>(args)...); }
    bool push(const T& v) { return try_enqueue(v); }
    bool push(T&& v) { return try_enqueue(std::move(v)); }

    bool try_dequeue(T& out) {
        std::size_t pos = dequeue_pos_.load(std::memory_order_relaxed);
        for (;;) {
            cell_t* cell = &buffer_[index(pos)];
            std::size_t seq = cell->seq.load(std::memory_order_acquire);
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
        for (;;) {
            if (try_dequeue(out)) {
                spaces_.release();
                return true;
            }

            if (closed_.load(std::memory_order_acquire)) {
                return try_dequeue(out);
            }
            items_.acquire();
        }
    }

    bool pop() {
        T dummy;
        return try_dequeue(dummy);
    }

    bool empty() const {
        return approx_size() == 0;
    }

    std::size_t capacity() const { return capacity_; }

    std::size_t approx_size() const {
        std::size_t enq = enqueue_pos_.load(std::memory_order_relaxed);
        std::size_t deq = dequeue_pos_.load(std::memory_order_relaxed);
        return enq >= deq ? (enq - deq) : 0;
    }

    void close() noexcept {
        closed_.store(true, std::memory_order_release);
        items_.release(consumer_count_);
    }

private:

    struct cell_t {
        std::atomic<std::size_t> seq;
        typename std::aligned_storage<sizeof(T), alignof(T)>::type storage;
        cell_t() : seq(0) {}
    };

    static constexpr bool is_power_of_two(std::size_t x) {
        return (x & (x - 1)) == 0;
    }

    static std::size_t round_up_to_power_of_two(std::size_t x) {
        if (is_power_of_two(x)) return x;
        --x;
        for (std::size_t i = 1; i < sizeof(std::size_t) * 8; i <<= 1)
            x |= x >> i;
        return x + 1;
    }

    static T* cell_value_ptr(cell_t* c) noexcept {
        return reinterpret_cast<T*>(&c->storage);
    }

    std::size_t index(std::size_t pos) const noexcept { return pos & mask_; }

    template <class... Args>
    bool emplace_impl(Args&&... args) {
        for (;;) {
            spaces_.acquire();
            std::size_t pos = enqueue_pos_.load(std::memory_order_relaxed);
            cell_t* cell = &buffer_[index(pos)];
            std::size_t seq = cell->seq.load(std::memory_order_acquire);
            intptr_t dif = static_cast<intptr_t>(seq) - static_cast<intptr_t>(pos);

            if (dif == 0) {
                new (cell_value_ptr(cell)) T(std::forward<Args>(args)...);
                cell->seq.store(pos + 1, std::memory_order_release);
                enqueue_pos_.store(pos + 1, std::memory_order_relaxed);
                items_.release();
                return true;
            }
            if (dif < 0) {
                //return false; // full
            }
        }
    }

    const std::size_t capacity_;
    const std::size_t mask_;
    const int consumer_count_;
    std::atomic<bool> closed_{ false };
    cell_t* const buffer_;
    alignas(64) std::atomic<std::size_t> enqueue_pos_;
    char _pad0[64 - sizeof(enqueue_pos_) % 64]{};
    alignas(64) std::atomic<std::size_t> dequeue_pos_;
    char _pad1[64 - sizeof(dequeue_pos_) % 64]{};
    std::counting_semaphore<> items_{ (ptrdiff_t)0 };
    std::counting_semaphore<> spaces_{ (ptrdiff_t)capacity_ };

};
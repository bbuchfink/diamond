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

#if defined(__APPLE__)
#include <dispatch/dispatch.h>
#elif defined(_WIN32)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <atomic>
#pragma comment(lib, "Synchronization.lib") // Required for WaitOnAddress
#elif defined(__linux__)
#include <unistd.h>
#include <sys/syscall.h>
#include <linux/futex.h>
#include <atomic>
#else
#error "Platform not supported"
#endif

template <ptrdiff_t least_max_value = 2147483647>
class CountingSemaphore {

#if defined(__APPLE__)
    dispatch_semaphore_t sem;
#else
    std::atomic<int32_t> count;
#endif

public:
    static constexpr ptrdiff_t max() noexcept { return least_max_value; }

    explicit CountingSemaphore(ptrdiff_t desired)
#if defined(__APPLE__)
        : sem(dispatch_semaphore_create(static_cast<intptr_t>(desired)))
#else
        : count(static_cast<int32_t>(desired))
#endif
    {
    }

#if defined(__APPLE__)
    ~CountingSemaphore() {
#if !__has_feature(objc_arc)
        dispatch_release(sem);
#endif
    }
#else
    ~CountingSemaphore() = default;
#endif

    CountingSemaphore(const CountingSemaphore&) = delete;
    CountingSemaphore& operator=(const CountingSemaphore&) = delete;

    void release(ptrdiff_t update = 1) {
#if defined(__APPLE__)        
        for (ptrdiff_t i = 0; i < update; ++i) {
            dispatch_semaphore_signal(sem);
        }
#else
        count.fetch_add(static_cast<int32_t>(update), std::memory_order_release);
        sys_wake(static_cast<int32_t>(update));
#endif
    }

    void acquire() {
#if defined(__APPLE__)
        dispatch_semaphore_wait(sem, DISPATCH_TIME_FOREVER);
#else
        int32_t old = count.load(std::memory_order_relaxed);
        while (true) {
            if (old > 0) {
                if (count.compare_exchange_weak(old, old - 1, std::memory_order_acquire, std::memory_order_relaxed)) {
                    return;
                }
            }
            else {
                sys_wait(old);
                old = count.load(std::memory_order_relaxed);
            }
        }
#endif
    }

private:

#if defined(_WIN32)

    void sys_wake(int32_t upd) {
        if (upd == 1) {
            WakeByAddressSingle(reinterpret_cast<void*>(&count));
        }
        else {
            WakeByAddressAll(reinterpret_cast<void*>(&count));
        }
    }

    void sys_wait(int32_t expected) {
        WaitOnAddress(reinterpret_cast<volatile void*>(&count), &expected, sizeof(int32_t), INFINITE);
    }
    
#elif defined(__linux__)

    void sys_wake(int32_t upd) {
        syscall(SYS_futex, reinterpret_cast<int32_t*>(&count), FUTEX_WAKE_PRIVATE, upd, nullptr, nullptr, 0);
    }

    void sys_wait(int32_t expected) {
        syscall(SYS_futex, reinterpret_cast<int32_t*>(&count), FUTEX_WAIT_PRIVATE, expected, nullptr, nullptr, 0);
    }

#endif

};
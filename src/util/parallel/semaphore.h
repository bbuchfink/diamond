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

#pragma once

// problems have been observed with std::counting_semaphore and clang as of version 18.1.8 on macOS
#if __cpp_lib_semaphore >= 201907L && (!defined(__clang__) || (__clang_major__ >= 21))
#include <semaphore>
#else
#warning "Old compiler missing semaphore support."

#include <mutex>
#include <condition_variable>

namespace std {

template <std::ptrdiff_t LeastMaxValue = -1>
class counting_semaphore {
private:
    std::ptrdiff_t counter;
    std::mutex mtx;
    std::condition_variable cv;

public:
    explicit counting_semaphore(std::ptrdiff_t desired = 0) : counter(desired) {}
    counting_semaphore(const counting_semaphore&) = delete;
    counting_semaphore& operator=(const counting_semaphore&) = delete;

    void release(std::ptrdiff_t update = 1) {
        {
            std::lock_guard<std::mutex> lock(mtx);
            counter += update;
        }
        if (update > 1) {
            cv.notify_all();
        } else {
            cv.notify_one();
        }
    }

    void acquire() {
        std::unique_lock<std::mutex> lock(mtx);
        cv.wait(lock, [this] { return counter > 0; });
        --counter;
    }
    
    static constexpr std::ptrdiff_t max() noexcept {
        return LeastMaxValue;
    }

};

using binary_semaphore = counting_semaphore<1>;

}

#endif
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
#if _MSC_VER || __GNUC__ >= 5 || __clang__
#define HAVE_SEMAPHORE __has_include(<semaphore>)
#endif

#if HAVE_SEMAPHORE
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
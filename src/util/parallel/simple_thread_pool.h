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

#include <vector>
#include <thread>
#include <atomic>
#include <mutex>
#include <functional>
#include <exception>

struct SimpleThreadPool {

    const std::atomic<bool>& stop() const {
        return stop_flag;
    }

    ~SimpleThreadPool() {
        stop_flag = true;
        for (auto& t : threads) {
            if (t.joinable()) {
                t.join();
            }
        }
    }

    template <typename Func>
    void spawn(Func&& func) {
        std::lock_guard<std::mutex> lock(data_mutex);

        threads.emplace_back([this, task = std::forward<Func>(func)]() {
            try {
                task(static_cast<const std::atomic<bool>&>(stop_flag));
            }
            catch (...) {
                std::lock_guard<std::mutex> lock(data_mutex);
                if (!first_exception) {
                    first_exception = std::current_exception();
                    stop_flag = true;
                }
            }
            });
    }

    void join_all() {
        for (auto& t : threads) {
            if (t.joinable()) {
                t.join();
            }
        }
        if (first_exception) {
            std::rethrow_exception(first_exception);
        }
    }

    void request_stop() {
        stop_flag = true;
    }

private:

    std::vector<std::thread> threads;
    std::atomic<bool> stop_flag{ false };
    std::exception_ptr first_exception = nullptr;
    std::mutex data_mutex;

};
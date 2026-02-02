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
#include <vector>
#include <thread>
#include <atomic>
#include <mutex>
#include <exception>
#include <functional>

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

    template <typename Func, typename... Args>
    std::thread::id spawn(Func&& func, Args&&... args) {
        std::lock_guard<std::mutex> lock(data_mutex);

        threads.emplace_back(
            [this,
            fn = std::forward<Func>(func),
            tup = std::make_tuple(std::forward<Args>(args)...)
            ]() mutable {
                try {
                    std::apply(
                        [&](auto&&... xs) {
                            std::invoke(fn,
                                static_cast<const std::atomic<bool>&>(stop_flag),
                                std::forward<decltype(xs)>(xs)...);
                        },
                        tup
                    );
                }
                catch (...) {
                    std::lock_guard<std::mutex> lock(data_mutex);
                    if (!first_exception) {
                        first_exception = std::current_exception();
                        stop_flag = true;
                    }
                }
            }
        );

        return threads.back().get_id();
    }

    template <typename T, typename Method, typename... Args>
    std::thread::id spawn_method(T* obj, Method&& method, Args&&... args) {
        if (!obj) {
            throw std::invalid_argument("spawn_method: obj must not be null");
        }
        return spawn(
            [obj, m = std::forward<Method>(method)](const std::atomic<bool>& stop, auto&&... xs) mutable {
                std::invoke(m, obj, stop, std::forward<decltype(xs)>(xs)...);
            },
            std::forward<Args>(args)...
        );
    }

    void join_all() {
        for (auto& t : threads) {
            if (t.joinable()) {
                t.join();
            }
        }
        threads.clear();
        if (first_exception) {
            std::rethrow_exception(first_exception);
        }
    }

    template<typename It>
    void join(It begin, It end) {
        for (It i = begin; i != end; ++i) {
            join(*i);
        }
        if (first_exception) {
            std::rethrow_exception(first_exception);
        }
    }

    void join(std::thread::id thread_id) {
        auto it = std::find_if(threads.begin(), threads.end(),
			[thread_id](const std::thread& t) { return t.get_id() == thread_id; });
        if(it == threads.end()) {
            throw std::runtime_error("Thread ID not found in thread pool.");
		}
        if (it->joinable()) {
            it->join();
		}
        threads.erase(it);
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
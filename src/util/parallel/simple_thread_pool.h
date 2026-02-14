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
#include <map>
#include <thread>
#include <atomic>
#include <mutex>
#include <exception>
#include <tuple>

namespace simple_thread_pool_detail {

template<std::size_t... Is>
struct index_sequence {};

template<std::size_t N, std::size_t... Is>
struct make_index_sequence_impl : make_index_sequence_impl<N-1, N-1, Is...> {};

template<std::size_t... Is>
struct make_index_sequence_impl<0, Is...> {
    typedef index_sequence<Is...> type;
};

template<std::size_t N>
using make_index_sequence = typename make_index_sequence_impl<N>::type;

template<typename Func, typename Tuple>
struct SpawnTask {
    Func fn;
    Tuple args;
    std::atomic<bool>* stop;
    std::mutex* mtx;
    std::exception_ptr* exc;

    SpawnTask(Func f, Tuple a, std::atomic<bool>* s, std::mutex* m, std::exception_ptr* e)
        : fn(std::move(f)), args(std::move(a)), stop(s), mtx(m), exc(e) {}

    void operator()() {
        try {
            invoke(make_index_sequence<std::tuple_size<Tuple>::value>{});
        }
        catch (...) {
            std::lock_guard<std::mutex> lock(*mtx);
            if (!*exc) {
                *exc = std::current_exception();
                *stop = true;
            }
        }
    }

    template<std::size_t... Is>
    void invoke(index_sequence<Is...>) {
        fn(static_cast<const std::atomic<bool>&>(*stop), std::get<Is>(args)...);
    }
};

template<typename T, typename Method>
struct MethodBinder {
    T* obj;
    Method method;

    MethodBinder(T* o, Method m) : obj(o), method(std::move(m)) {}

    template<typename... Xs>
    void operator()(const std::atomic<bool>& stop, Xs&&... xs) {
        (obj->*method)(stop, std::forward<Xs>(xs)...);
    }
};

}

struct SimpleThreadPool {

    const std::atomic<bool>& stop() const {
        return stop_flag;
    }

    ~SimpleThreadPool() {
        stop_flag = true;
        for (auto& t : threads) {
            if (t.second.joinable()) {
                t.second.join();
            }
        }
    }

    template <typename Func, typename... Args>
    std::thread::id spawn(Func&& func, Args&&... args) {
        std::lock_guard<std::mutex> lock(data_mutex);

        using FuncD = typename std::decay<Func>::type;
        using TupleD = std::tuple<typename std::decay<Args>::type...>;
        using Task = simple_thread_pool_detail::SpawnTask<FuncD, TupleD>;

        std::thread t(Task(
            std::forward<Func>(func),
            std::make_tuple(std::forward<Args>(args)...),
            &stop_flag, &data_mutex, &first_exception
        ));

		std::thread::id id = t.get_id();
        threads[id] = std::move(t);
        return id;
    }

    template <typename T, typename Method, typename... Args>
    std::thread::id spawn_method(T* obj, Method&& method, Args&&... args) {
        if (!obj) {
            throw std::invalid_argument("spawn_method: obj must not be null");
        }
        using MethodD = typename std::decay<Method>::type;
        return spawn(
            simple_thread_pool_detail::MethodBinder<T, MethodD>(obj, std::forward<Method>(method)),
            std::forward<Args>(args)...
        );
    }

    void join_all() {
        for (auto& t : threads) {
            if (t.second.joinable()) {
                t.second.join();
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
        std::thread t;
        {
            std::lock_guard<std::mutex> lock(data_mutex);
            auto it = threads.find(thread_id);
            if (it == threads.end()) {
                throw std::runtime_error("Thread ID not found in thread pool.");
            }
            t = std::move(it->second);
            threads.erase(it);
        }        
        if (t.joinable()) {
            t.join();
		}
	}

    void request_stop() {
        stop_flag = true;
    }

private:

    std::map<std::thread::id, std::thread> threads;
    std::atomic<bool> stop_flag{ false };
    std::exception_ptr first_exception = nullptr;
    std::mutex data_mutex;

};
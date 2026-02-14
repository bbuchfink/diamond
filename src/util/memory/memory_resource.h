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
#if (__GNUC__ >= 5 && __has_include(<memory_resource>)) || (defined(__cpp_lib_memory_resource) && __cpp_lib_memory_resource >= 201603) || _MSC_VER
#include <memory_resource>
#else
#warning "Old compiler missing memory_resource support."
#include <list>
#include <vector>

namespace std { namespace pmr {

struct monotonic_buffer_resource {};

template<typename T>
struct list : public std::list<T> {
    list(monotonic_buffer_resource*):
        std::list<T>()
    { }
};

template<typename T>
struct vector : public std::vector<T> {
    vector(monotonic_buffer_resource*) :
        std::vector<T>()
    {
    }
};

}}

#endif
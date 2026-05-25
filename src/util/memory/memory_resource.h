/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#if __GNUC__ >= 5 || _MSC_VER || __clang__
#if __has_include(<memory_resource>)
#define HAVE_MEMORY_RESOURCE
#endif
#endif

#if defined(__cpp_lib_memory_resource) && __cpp_lib_memory_resource >= 201603
#define HAVE_MEMORY_RESOURCE
#endif

#ifdef HAVE_MEMORY_RESOURCE
#include <memory_resource>
#else
#warning "Old compiler missing memory_resource support."
#include <list>
#include <vector>
#include <string>
#include <unordered_map>

namespace std { namespace pmr {

struct memory_resource {};
struct monotonic_buffer_resource : public memory_resource {};
struct unsynchronized_pool_resource : public memory_resource {};

template<typename T>
struct list : public std::list<T> {
    list(memory_resource*) :
        std::list<T>()
    {
    }
    list(monotonic_buffer_resource*):
        std::list<T>()
    { }
};

template<typename T>
struct vector : public std::vector<T> {
    vector(size_t n, memory_resource*) :
        std::vector<T>(n)
    {
    }
    vector(memory_resource*) :
        std::vector<T>()
    {
    }
    vector(monotonic_buffer_resource*) :
        std::vector<T>()
    {
    }
};

struct string : public std::string {
    string():
        std::string()
    { }
    string(memory_resource*) :
        std::string()
    {
    }
    string(monotonic_buffer_resource*) :
        std::string()
    {
    }
    string(std::pmr::string&& s) :
        std::string(std::move(s))
    {
	}
    string(const std::pmr::string& s) :
        std::string(s)
    {
    }
    string(const std::pmr::string& s, memory_resource*) :
        std::string(s)
    {
    }
    string(const std::pmr::string& s, monotonic_buffer_resource*) :
        std::string(s)
    {
    }
    string& operator=(const std::string& s) {
        std::string::operator=(s);
        return *this;
    }
};

template<typename KT, typename VT>
struct unordered_map : public std::unordered_map<KT, VT> {
    unordered_map(memory_resource*) :
        std::unordered_map<KT, VT>()
    {
    }
    unordered_map(monotonic_buffer_resource*) :
        std::unordered_map<KT, VT>()
    {
    }
};

}}

#endif
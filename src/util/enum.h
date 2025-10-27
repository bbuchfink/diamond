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
#include <map>
#include <string>
#include <stdexcept>

template<typename T>
struct FieldValue{
	FieldValue(const T v, const bool primary = true):
		v(v),
		primary(primary)
	{}
	const T v;
	const bool primary;
};

template<typename T>
using EMap = std::map<T, std::string>;
template<typename T>
using SEMap = std::map<std::string, FieldValue<T>>;

template<typename T> struct EnumTraits {};

template<typename T> std::string to_string(T v) {
	auto it = EnumTraits<T>::to_string.find(v);
	if (it == EnumTraits<T>::to_string.end())
		throw std::runtime_error("Invalid conversion from enum to string.");
	return it->second;
}

template<typename T>
T from_string(const std::string& s) {
	auto it = EnumTraits<T>::from_string.find(s);
	if (it == EnumTraits<T>::from_string.end()) {
		std::string p;
		size_t n = 0;
		for (const auto& i : EnumTraits<T>::from_string) {
			if (!i.second.primary)
				continue;
			if (n++)
				p += ", ";
			p += i.first;
		}
		throw std::runtime_error("Invalid value for string field: " + s + ". Permitted values: " + p);
	}
	return it->second.v;
}

#ifndef DEFINE_ENUM_FLAG_OPERATORS
#define DEFINE_ENUM_FLAG_OPERATORS(T) constexpr inline T operator~ (T a) { return static_cast<T>(~static_cast<std::underlying_type<T>::type>(a)); } \
constexpr inline T operator| (T a, T b) { return static_cast<T>(static_cast<std::underlying_type<T>::type>(a) | static_cast<std::underlying_type<T>::type>(b)); } \
constexpr inline T operator& (T a, T b) { return static_cast<T>(static_cast<std::underlying_type<T>::type>(a) & static_cast<std::underlying_type<T>::type>(b)); } \
constexpr inline T operator^ (T a, T b) { return static_cast<T>(static_cast<std::underlying_type<T>::type>(a) ^ static_cast<std::underlying_type<T>::type>(b)); } \
inline T& operator|= (T& a, T b) { return reinterpret_cast<T&>(reinterpret_cast<std::underlying_type<T>::type&>(a) |= static_cast<std::underlying_type<T>::type>(b)); } \
inline T& operator&= (T& a, T b) { return reinterpret_cast<T&>(reinterpret_cast<std::underlying_type<T>::type&>(a) &= static_cast<std::underlying_type<T>::type>(b)); } \
inline T& operator^= (T& a, T b) { return reinterpret_cast<T&>(reinterpret_cast<std::underlying_type<T>::type&>(a) ^= static_cast<std::underlying_type<T>::type>(b)); }
#endif

template<typename T>
bool flag_all(T a, T b) {
	typedef typename std::underlying_type<T>::type U;
	U c = static_cast<U>(b);
	return (static_cast<U>(a) & c) == c;
}

template<typename T>
bool flag_any(T a, T b) {
	typedef typename std::underlying_type<T>::type U;
	U c = static_cast<U>(b);
	return (static_cast<U>(a) & c) != U(0);
}

template<typename T>
bool flag_only(T a, T b) {
	typedef typename std::underlying_type<T>::type U;
	U c = static_cast<U>(b);
	return (static_cast<U>(a) & ~c) == U(0);
}
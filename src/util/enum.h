/****
DIAMOND protein aligner
Copyright (C) 2021 Max Planck Society for the Advancement of Science e.V.

Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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

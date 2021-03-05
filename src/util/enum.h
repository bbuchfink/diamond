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

template<typename _t>
using EMap = std::map<_t, std::string>;

template<typename _t> struct EnumTraits {};

template<typename _t> std::string to_string(_t v) {
	auto it = EnumTraits<_t>::to_string.find(v);
	if (it == EnumTraits<_t>::to_string.end())
		throw std::runtime_error("Invalid conversion from enum to string.");
	return it->second;
}

#define DEF_ENUM_FLAG_OPERATORS(T) constexpr inline T operator~ (T a) { return static_cast<T>(~static_cast<std::underlying_type<T>::type>(a)); } \
constexpr inline T operator| (T a, T b) { return static_cast<T>(static_cast<std::underlying_type<T>::type>(a) | static_cast<std::underlying_type<T>::type>(b)); } \
constexpr inline T operator& (T a, T b) { return static_cast<T>(static_cast<std::underlying_type<T>::type>(a) & static_cast<std::underlying_type<T>::type>(b)); } \
constexpr inline T operator^ (T a, T b) { return static_cast<T>(static_cast<std::underlying_type<T>::type>(a) ^ static_cast<std::underlying_type<T>::type>(b)); } \
inline T& operator|= (T& a, T b) { return reinterpret_cast<T&>(reinterpret_cast<std::underlying_type<T>::type&>(a) |= static_cast<std::underlying_type<T>::type>(b)); } \
inline T& operator&= (T& a, T b) { return reinterpret_cast<T&>(reinterpret_cast<std::underlying_type<T>::type&>(a) &= static_cast<std::underlying_type<T>::type>(b)); } \
inline T& operator^= (T& a, T b) { return reinterpret_cast<T&>(reinterpret_cast<std::underlying_type<T>::type&>(a) ^= static_cast<std::underlying_type<T>::type>(b)); }

template<typename _t>
bool flag_get(_t a, _t b) {
	typedef typename std::underlying_type<_t>::type T;
	T c = static_cast<T>(b);
	return static_cast<T>(a) & c == c;
}
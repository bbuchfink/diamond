/****
DIAMOND protein aligner
Copyright (C) 2024 Benjamin Buchfink

Code developed by Benjamin Buchfink <buchfink@gmail.com>

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

#include "util/algo/varint.h"
#include "util/io/serializer.h"

static inline void write_varint(Serializer& s, int32_t x) {
	char buf[5];
	char* end = write_varuint32(x, buf);
	s.write_raw(buf, end - buf);
}

template<typename It>
static inline void serialize_varint(Serializer& s, It begin, It end) {
	for (It i = begin; i != end; ++i)
		write_varint(s, *i);
}

template<typename It>
static inline void serialize(Serializer& s, It begin, It end) {
	for (It i = begin; i != end; ++i)
		s << *i;
}

static inline void serialize(Serializer& s, const std::set<int32_t>& v) {
	write_varint(s, (int32_t)v.size());
	serialize_varint(s, v.cbegin(), v.cend());
}

static inline void serialize(Serializer& s, const std::vector<int32_t>& v) {
	s.write((uint32_t)v.size());
	serialize(s, v.cbegin(), v.cend());
}

static inline void serialize(Serializer& s, const std::vector<std::string>& v) {
	s.write((uint32_t)v.size());
	serialize(s, v.cbegin(), v.cend());
}

template<typename T1, typename T2>
void serialize(Serializer& s, const std::pair<T1, T2>& p) {
	s << p.first << p.second;
}

static inline void deserialize(Deserializer& d, std::vector<std::string>& out) {
	uint32_t n;
	d >> n;
	out.clear();
	out.reserve(n);
	std::string s;
	for (uint32_t i = 0; i < n; ++i) {
		d >> s;
		out.push_back(std::move(s));
	}
}

static inline void deserialize(Deserializer& d, std::vector<std::int32_t>& out) {
	uint32_t n;
	d >> n;
	out.clear();
	out.reserve(n);
	int32_t x;
	for (uint32_t i = 0; i < n; ++i) {
		d >> x;
		out.push_back(x);
	}
}

template<typename T1, typename T2>
void deserialize(Deserializer& s, std::pair<T1, T2>& out) {
	s >> out.first >> out.second;
}
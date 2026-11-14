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
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>
#include "util/parallel/filestack.h"

// Parses one line of a list file: a path followed by up to max_fields unsigned integer fields.
// The files written by the clustering workflow separate fields by tabs, so the path may contain
// spaces. Lines without a tab (e.g. user supplied volume lists) are split on whitespace, the path
// being everything before the trailing integer fields. Returns false for blank lines.
inline bool parse_list_line(const std::string& line, std::string& path, std::vector<uint64_t>& fields, size_t max_fields) {
	static const char* const WS = " \t\r\n";
	fields.clear();
	const size_t b = line.find_first_not_of(WS);
	if (b == std::string::npos)
		return false;
	const size_t e = line.find_last_not_of(WS) + 1;
	std::string s = line.substr(b, e - b);
	auto is_uint = [](const std::string& t) {
		return !t.empty() && t.size() <= 20 && std::all_of(t.begin(), t.end(), [](char c) { return c >= '0' && c <= '9'; });
	};
	const size_t tab = s.find('\t');
	if (tab != std::string::npos) {
		path = s.substr(0, tab);
		size_t i = tab + 1;
		while (i <= s.size()) {
			size_t j = s.find('\t', i);
			if (j == std::string::npos)
				j = s.size();
			std::string t = s.substr(i, j - i);
			const size_t tb = t.find_first_not_of(WS);
			t = tb == std::string::npos ? std::string() : t.substr(tb, t.find_last_not_of(WS) + 1 - tb);
			if (!t.empty()) {
				if (!is_uint(t))
					throw std::runtime_error("Format error in list file line: " + line);
				if (fields.size() < max_fields)
					fields.push_back(std::stoull(t));
			}
			i = j + 1;
		}
	}
	else {
		while (fields.size() < max_fields) {
			const size_t sp = s.find_last_of(WS);
			if (sp == std::string::npos)
				break;
			const std::string t = s.substr(sp + 1);
			if (!is_uint(t))
				break;
			fields.push_back(std::stoull(t));
			s.resize(s.find_last_not_of(WS, sp) + 1);
		}
		std::reverse(fields.begin(), fields.end());
		path = s;
	}
	return true;
}

struct Bucket {
	static constexpr uint64_t NIL = std::numeric_limits<uint64_t>::max();
	Bucket() = default;
	Bucket(const std::string& path, uint64_t records, uint64_t key_begin = NIL, uint64_t key_end = NIL) : path(path), records_(records), key_begin_(key_begin), key_end_(key_end) {}
	std::string path;
	std::string containing_directory() const {
		return ::containing_directory(path);
	}
	uint64_t key_begin() const {
		if (key_begin_ == NIL)
			throw std::runtime_error("Key range not set for bucket: " + path);
		return key_begin_;
	}
	uint64_t key_end() const {
		if (key_end_ == NIL)
			throw std::runtime_error("Key range not set for bucket: " + path);
		return key_end_;
	}
	uint64_t records() const {
		if (records_ == NIL)
			throw std::runtime_error("Record count not set for bucket: " + path);
		return records_;
	}
private:
	uint64_t records_ = NIL, key_begin_ = NIL, key_end_ = NIL;
	friend struct RadixedTable;
};

struct RadixedTable : public std::vector<Bucket> {

	RadixedTable(int shift):
		shift(shift)
	{
	}

	RadixedTable(const std::string& file_name, int shift) : shift(shift) {
		std::ifstream file(file_name);
		std::string line, path;
		std::vector<uint64_t> fields;
		while (std::getline(file, line)) {
			if (!parse_list_line(line, path, fields, 1))
				continue;
			if (fields.empty())
				throw std::runtime_error("Format error in RadixedTable");
			emplace_back(path, fields.front());
		}
	}

	uint64_t max_buckets(uint64_t mem_limit, size_t record_size) const {
		std::vector<uint64_t> counts;
		counts.reserve(size());
		for (const Bucket& b : *this)
			counts.push_back(b.records());
		std::sort(counts.begin(), counts.end(), std::greater<uint64_t>());
		uint64_t sum = 0;
		for (uint64_t i = 0; i < counts.size(); ++i) {
			sum += counts[i] * record_size;
			if (sum >= mem_limit)
				return i > 0 ? i : 1;
		}
		return counts.size();
	}

	void append(FileStack& out) const {
		std::ostringstream ss;
		for (const Bucket& b : *this) {
			ss << b.path << '\t' << b.records_ << std::endl;
		}
		out.push(ss.str());
	}

	const int shift;

};
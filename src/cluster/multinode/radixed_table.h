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

#include <vector>
#include "util/parallel/filestack.h"

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
		std::string path;
		uint64_t records;
		for (;;) {
			if (!(file >> path))
				return;
			file >> records;
			if (!file)
				throw std::runtime_error("Format error in RadixedTable");
			emplace_back(path, records);
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
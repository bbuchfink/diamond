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
#include <fstream>
#include "util/system/system.h"

namespace External {

struct Volume {
	Volume() :
		path(),
		oid_begin(0),
		oid_end(0),
		record_count(0)
	{
	}
	Volume(const std::string& path, OId oid_begin, OId oid_end, OId record_count) :
		path(path),
		oid_begin(oid_begin),
		oid_end(oid_end),
		record_count(record_count)
	{
	}
	std::string path;
	OId oid_begin, oid_end, record_count;
	bool operator<(size_t oid) const {
		return oid_end <= oid;
	}
	OId oid_range() const {
		return oid_end - oid_begin;
	}
	bool operator<(const Volume& v) const {
		return oid_begin < v.oid_begin;
	}
	friend std::istream& operator>>(std::istream& str, Volume& v) {
		std::string line;
		v.oid_begin = v.oid_end = 0;
		if (!std::getline(str, line)) return str;
		std::istringstream row(line);
		row >> v.path;
		if (!row)
			throw std::runtime_error("Format error in VolumedFile");
		row >> v.record_count;
		if (!row)
			throw std::runtime_error("Format error in VolumedFile");
		row >> v.oid_begin >> v.oid_end;
		return str;
	}
};

struct Bucket {
	static constexpr uint64_t NIL = std::numeric_limits<uint64_t>::max();
	Bucket() = default;
	Bucket(const std::string& path, uint64_t records) : path(path), records_(records) {}
	std::string path;
	std::string containing_directory() const {
		return ::containing_directory(path);
	}
	uint64_t records() const {
		if (records_ == NIL)
			throw std::runtime_error("Record count not set for bucket: " + path);
		return records_;
	}
private:
	uint64_t records_ = NIL;
	friend struct RadixedTable;
};

struct VolumedFile : public std::vector<Volume> {
	VolumedFile(const Bucket& bucket):
		VolumedFile(bucket.path)
	{ }
	VolumedFile(const std::string& file_name) :
		list_file_(file_name),
		records_(0),
		max_oid_(0)
	{
		std::ifstream volume_file(file_name);
		if (!volume_file)
			throw std::runtime_error("Error opening file " + file_name);
		int64_t oid = 0;
		Volume v;
		while (volume_file) {
			volume_file >> v;
			if (!volume_file)
				break;
			if (v.oid_begin == 0 && v.oid_end == 0) {
				v.oid_begin = oid;
				v.oid_end = oid + v.record_count;
			}
			push_back(v);
			oid += v.record_count;
			records_ += v.record_count;
			max_oid_ = std::max(max_oid_, v.oid_end - 1);
		}
		std::sort(begin(), end());
	}
	OId sparse_records() const {
		return records_;
	}
	OId max_oid() const {
		return max_oid_;
	}
	std::pair<std::vector<Volume>::const_iterator, std::vector<Volume>::const_iterator> find(OId oid_begin, OId oid_end) const {
		auto it = std::lower_bound(begin(), end(), oid_begin);
		if (it == end())
			throw std::runtime_error("OID out of bounds");
		auto end = it + 1;
		while (end < this->end() && end->oid_begin < oid_end)
			++end;
		return { it,end };
	}
	void remove() const {
		for (const Volume& v : *this)
			::remove(v.path.c_str());
		::remove(list_file_.c_str());
		rmdir(containing_directory(list_file_).c_str());
	}
private:
	const std::string list_file_;
	OId records_, max_oid_;
};

struct RadixedTable : public std::vector<Bucket> {

	RadixedTable() = default;

	RadixedTable(const std::string& file_name) {
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
		for(const Bucket& b : *this) {
			ss << b.path << '\t' << b.records() << std::endl;
		}
		out.push(ss.str());
	}

};

}
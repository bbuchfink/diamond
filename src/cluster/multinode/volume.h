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
#include <algorithm>
#include <sstream>
#include <memory>
#include "util/system/system.h"
#include "radixed_table.h"
#include "data/sequence_file.h"
#include "data/fasta/fasta_file.h"

struct Volume {
	Volume() :
		path(),
		oid_begin(0),
		oid_end(0),
		record_count(std::numeric_limits<OId>::max())
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
		row >> v.oid_begin >> v.oid_end;
		return str;
	}
};

struct VolumedFile : public std::vector<Volume> {
	VolumedFile(const Bucket& bucket) :
		VolumedFile(bucket.path)
	{
	}
	VolumedFile(const std::string& file_name, uint64_t letter_count = 0) :
		list_file_(file_name),
		records_(std::numeric_limits<OId>::max()),
		max_oid_(0),
		letter_count_(letter_count)
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
			if (v.record_count == std::numeric_limits<OId>::max()) {
				if (records_ != std::numeric_limits<OId>::max())
					throw std::runtime_error("Inconsistent record count");
			}
			else {
				if (records_ == std::numeric_limits<OId>::max())
					records_ = 0;
				records_ += v.record_count;
			}
			if (v.oid_end > 0)
				max_oid_ = std::max(max_oid_, v.oid_end - 1);
		}
		std::sort(begin(), end());
		if (empty())
			records_ = 0;
	}
	OId sparse_records() const {
		if (records_ == std::numeric_limits<OId>::max())
			throw std::runtime_error("Record count not set");
		return records_;
	}
	OId max_oid() const {
		if (max_oid_ == 0)
			throw std::runtime_error("max_oid not set");
		return max_oid_;
	}
	void set_max_oid(OId max_oid) {
		max_oid_ = max_oid;
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
	std::vector<Volume>::const_iterator find(OId oid) const {
		auto it = std::lower_bound(begin(), end(), oid);
		if (it == end() || oid < it->oid_begin || oid >= it->oid_end)
			throw std::runtime_error("OID out of bounds");
		return it;
	}
	void remove(bool dir = true, bool files = true) const {
		if (files)
			for (const Volume& v : *this)
				remove_tmp_file(v.path);
		remove_tmp_file(list_file_);
		if (dir)
			rmdir(containing_directory(list_file_).c_str());
	}
	void set_letter_count(uint64_t count) {
		letter_count_ = count;
	}
	uint64_t letter_count() const {
		if (letter_count_ == 0)
			throw std::runtime_error("Letter count not set");
		return letter_count_;
	}
	const std::string& list_file() const {
		return list_file_;
	}
private:
	const std::string list_file_;
	OId records_, max_oid_;
	uint64_t letter_count_;
};

inline Block* load_seqs(const VolumedFile& volumes, OId oid_begin, OId oid_end, const std::string& index_dir, SequenceFile::Flags flags = SequenceFile::Flags::ALL) {
	std::vector<Volume>::const_iterator vol_begin, vol_end;
	std::tie(vol_begin, vol_end) = volumes.find(oid_begin, oid_end);
	Block* combined = nullptr;
	for (auto v = vol_begin; v != vol_end; ++v) {
		const OId local_begin = std::max(oid_begin, v->oid_begin) - v->oid_begin;
		const OId local_end = std::min(oid_end, v->oid_end) - v->oid_begin;
		const OId count = local_end - local_begin;
		//std::unique_ptr<SequenceFile> file(SequenceFile::auto_create({ v->path }, flags));
		std::unique_ptr<SequenceFile> file(new FastaFile({ v->path }, flags, amino_acid_traits, index_dir + std::to_string(v - volumes.begin())));
		Block* vol_block;
		file->set_seqinfo_ptr(local_begin);
		vol_block = file->load_seqs(INT64_MAX, count);
		vol_block->offset_oids(v->oid_begin);
		if (!combined) {
			combined = vol_block;
		}
		else {
			combined->append(*vol_block, true);
			delete vol_block;
		}
		file->close();
	}
	if (!combined)
		combined = new Block();
	return combined;
}
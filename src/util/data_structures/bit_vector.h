/****
DIAMOND protein aligner
Copyright (C) 2020 Max Planck Society for the Advancement of Science e.V.

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

#include <mutex>
#include <vector>
#include <stdint.h>
#include <algorithm>
#include "../intrin.h"

struct BitVector {

	BitVector() :
		size_(0)
	{}

	BitVector(size_t n) :
		data_((n + 63) / 64, 0),
		size_(n)
	{
	}

	void set(size_t i) {
		data_[i >> 6] |= uint64_t(1) << (i & 63);
	}

	bool get(size_t i) const {
		return data_[i >> 6] & (uint64_t(1) << (i & 63));
	}

	BitVector& operator|=(const BitVector& v) {
		for (size_t i = 0; i < data_.size(); ++i)
			data_[i] |= v.data_[i];
		return *this;
	}

	void reset() {
		std::fill(data_.begin(), data_.end(), 0);
	}

	size_t one_count() const {
		size_t n = 0;
		for (uint64_t x : data_)
			n += popcount64(x);
		return n;
	}

	bool empty() const {
		return data_.empty();
	}

	int64_t size() const {
		return size_;
	}

	std::vector<int64_t> negative_list() const {
		std::vector<int64_t> v;
		for (int64_t i = 0; i < size_; ++i)
			if (!get(i))
				v.push_back(i);
		return v;
	}

private:

	std::vector<uint64_t> data_;
	int64_t size_;

};
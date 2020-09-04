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

struct BitVector {

	BitVector(size_t n) :
		data_((n + 63) / 64, 0)
	{
	}

	void set(size_t i) {
		data_[i >> 6] |= uint64_t(1) << (i & 63);
	}

	bool get(size_t i) {
		return data_[i >> 6] & (uint64_t(1) << (i & 63));
	}

	BitVector& operator|=(const BitVector& v) {
		std::lock_guard<std::mutex> lock(mtx_);
		for (size_t i = 0; i < data_.size(); ++i)
			data_[i] |= v.data_[i];
		return *this;
	}

private:

	std::vector<uint64_t> data_;
	std::mutex mtx_;

};
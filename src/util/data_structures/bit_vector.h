/****
DIAMOND protein aligner
Copyright (C) 2020 Max Planck Society for the Advancement of Science e.V.

Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>
Arm NEON port contributed by Martin Larralde <martin.larralde@embl.de>

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
#include "../simd.h"

struct BitVector {

	BitVector() {}

	BitVector(size_t n) :
		data_((n + 63) / 64, 0)
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
		size_t i = 0;
		size_t n = 0;
		const uint64_t* data = data_.data();

#ifdef __ARM_NEON
		uint16x8_t acc;
		for (; i + 2 < data_.size(); ++i) {
			const uint8x16_t block = vreinterpretq_u8_u64(vld1q_u64(&data[i]));
			const uint8x16_t count = vcntq_u8(block);
			acc = vpadalq_u8(acc, count);
			if ((i % std::numeric_limits<uint16_t>::max()) == std::numeric_limits<uint16_t>::max()) {
				n += ::SIMD::vhsumq_u16(acc);
				acc = veorq_u16(acc, acc);
			}
		}
		n += ::SIMD::vhsumq_u16(acc);
#endif
		for (; i < data_.size(); ++i)
			n += popcount64(data[i]);
		return n;
	}

	bool empty() const {
		return data_.empty();
	}

private:

	std::vector<uint64_t> data_;

};

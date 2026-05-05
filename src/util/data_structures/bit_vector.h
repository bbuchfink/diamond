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

#include <atomic>
#include <vector>
#include <stdint.h>
#include "../intrin.h"
#include "../simd.h"

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

	void atomic_set(size_t i) {
		const uint64_t mask = uint64_t(1) << (i & 63);
#if defined(__cpp_lib_atomic_ref)
		std::atomic_ref<uint64_t> word(data_[i >> 6]);
		word.fetch_or(mask, std::memory_order_seq_cst);
#elif defined(_MSC_VER)
		_InterlockedOr64(reinterpret_cast<volatile long long*>(&data_[i >> 6]), static_cast<long long>(mask));
#else
		__atomic_fetch_or(&data_[i >> 6], mask, __ATOMIC_SEQ_CST);
#endif
	}

	bool get(size_t i) const {
		return data_[i >> 6] & (uint64_t(1) << (i & 63));
	}

	bool atomic_load(size_t i) const {
#if defined(__cpp_lib_atomic_ref)
		const std::atomic_ref<uint64_t> atomic_word(const_cast<uint64_t&>(data_[i >> 6]));
		const uint64_t word = atomic_word.load(std::memory_order_seq_cst);
#elif defined(_MSC_VER)
		const uint64_t word = static_cast<uint64_t>(_InterlockedCompareExchange64(reinterpret_cast<volatile long long*>(&const_cast<uint64_t&>(data_[i >> 6])), 0, 0));
#else
		const uint64_t word = __atomic_load_n(&data_[i >> 6], __ATOMIC_SEQ_CST);
#endif
		return word & (uint64_t(1) << (i & 63));
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
			/* Each loop increments each lane by up to 8, so the accumulator (which has 16-bit lanes) 
			 * may overflow after 2^16 / 8 = 2^13 iterations.
			 */
			if ((i % (1 << 13)) == 0) {
				n += vhsumq_u16(acc);
				acc = veorq_u16(acc, acc);
			}
		}
		n += vhsumq_u16(acc);
#else
		for (; i < data_.size(); ++i)
			n += popcount64(data[i]);
#endif
		return n;
	}

	bool empty() const {
		return data_.empty();
	}

	uint64_t size() const {
		return size_;
	}

	std::vector<uint64_t> negative_list() const {
		std::vector<uint64_t> v;
		for (uint64_t i = 0; i < size_; ++i)
			if (!get(i))
				v.push_back(i);
		return v;
	}

private:

	std::vector<uint64_t> data_;
	uint64_t size_;

};
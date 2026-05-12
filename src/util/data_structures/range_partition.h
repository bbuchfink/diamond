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
#include <algorithm>
#include <utility>
#include <assert.h>
#include <limits.h>
#include <array>
#include "util/algo/sort.h"

template<int N, typename T>
struct RangePartition {

	RangePartition(const int *begin, int count, int end) :
		count_(1)
	{
		assert(count > 0);
		std::array<std::pair<int, int>, N> b;
		for (int i = 0; i < count; ++i)
			b[i] = { begin[i], i };
		insertion_sort(b.begin(), b.begin() + count);
		std::fill(mask_[0], mask_[0] + N, std::numeric_limits<T>::min());
		mask_[0][b[0].second] = 0;
#ifdef DP_STAT
		bit_mask_[0] = 1llu << b[0].second;
#endif
		begin_[0] = b[0].first;
		for (int i = 1; i < count; ++i) {
			if (begin_[count_ - 1] < b[i].first) {
				begin_[count_] = b[i].first;
				std::copy(mask_[count_ - 1], mask_[count_ - 1] + N, mask_[count_]);
				mask_[count_][b[i].second] = 0;
#ifdef DP_STAT
				bit_mask_[count_] = bit_mask_[count_ - 1] | (1llu << b[i].second);
#endif
				++count_;
			}
			else {
				mask_[count_ - 1][b[i].second] = 0;
#ifdef DP_STAT
				bit_mask_[count_ - 1] |= 1llu << b[i].second;
#endif
			}
		}
		begin_[count_] = end;
	}

	int begin(int i) const {
		return begin_[i];
	}

	int end(int i) const {
		return begin_[i + 1];
	}

	int count() const {
		return count_;
	}

	const T* mask(int i) const {
		return mask_[i];
	}

#ifdef DP_STAT
	uint64_t bit_mask(int i) const {
		return bit_mask_[i];
	}
#endif

private:

	int begin_[N + 1];
	T mask_[N][N];
#ifdef DP_STAT
	uint64_t bit_mask_[N];
#endif
	int count_;
	
};
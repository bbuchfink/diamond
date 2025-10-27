/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

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
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
#include <stdint.h>
#include <mutex>

struct Statistics
{

	using StatType = int64_t;

	enum value {
		SEED_HITS, TENTATIVE_MATCHES0, TENTATIVE_MATCHES1, TENTATIVE_MATCHES2, TENTATIVE_MATCHES3, TENTATIVE_MATCHES4, TENTATIVE_MATCHESX, MATCHES, ALIGNED, GAPPED, DUPLICATES,
		GAPPED_HITS, QUERY_SEEDS, QUERY_SEEDS_HIT, REF_SEEDS, REF_SEEDS_HIT, QUERY_SIZE, REF_SIZE, OUT_HITS, OUT_MATCHES, COLLISION_LOOKUPS, QCOV, BIAS_ERRORS, SCORE_TOTAL, ALIGNED_QLEN, PAIRWISE, HIGH_SIM,
		SEARCH_TEMP_SPACE, SECONDARY_HITS, ERASED_HITS, SQUARED_ERROR, CELLS, TARGET_HITS0, TARGET_HITS1, TARGET_HITS2, TARGET_HITS3, TARGET_HITS3_CBS, TARGET_HITS4, TARGET_HITS5, TARGET_HITS6, TIME_GREEDY_EXT, LOW_COMPLEXITY_SEEDS,
		SWIPE_REALIGN, EXT8, EXT16, EXT32, GAPPED_FILTER_TARGETS, GAPPED_FILTER_HITS1, GAPPED_FILTER_HITS2, GROSS_DP_CELLS, NET_DP_CELLS, TIME_TARGET_SORT, TIME_SW, TIME_EXT, TIME_GAPPED_FILTER,
		TIME_LOAD_HIT_TARGETS, TIME_CHAINING, TIME_LOAD_SEED_HITS, TIME_SORT_SEED_HITS, TIME_SORT_TARGETS_BY_SCORE, TIME_TARGET_PARALLEL, TIME_TRACEBACK_SW, TIME_TRACEBACK, HARD_QUERIES, TIME_MATRIX_ADJUST,
		MATRIX_ADJUST_COUNT, COMP_BASED_STATS_COUNT, FAILED_COMP_BASED_STATS, MASKED_LAZY, SWIPE_TASKS_TOTAL, SWIPE_TASKS_ASYNC, TRIVIAL_ALN, TIME_EXT_32, EXT_OVERFLOW_8, EXT_WASTED_16, DP_CELLS_8, DP_CELLS_16, DP_CELLS_32, TIME_PROFILE, TIME_ANCHORED_SWIPE,
		TIME_ANCHORED_SWIPE_ALLOC, TIME_ANCHORED_SWIPE_SORT, TIME_ANCHORED_SWIPE_ADD, TIME_ANCHORED_SWIPE_OUTPUT, TIME_PROFILE_GENERATION, EXTENSIONS_RECOMPUTE,
		TIME_SEARCH, SEEDS_HIT, COUNT
	};

	Statistics()
	{
		reset();
	}

	void reset() {
		std::fill(data_, data_ + COUNT, (StatType)0);
	}

	Statistics& operator+=(const Statistics &rhs)
	{
		mtx_.lock();
		for(unsigned i=0;i<COUNT;++i)
			data_[i] += rhs.data_[i];
		mtx_.unlock();
		return *this;
	}

	void inc(const value v, StatType n = 1lu)
	{ data_[v] += n; }

	void max(const value v, StatType n)
	{
		data_[v] = std::max(data_[v], n);
	}

	StatType get(const value v) const
	{ return data_[v]; }

	void print() const;

	StatType data_[COUNT];
	std::mutex mtx_;

};

extern Statistics statistics;
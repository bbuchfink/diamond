/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen
						
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
#include <algorithm>
#include <stdint.h>
#include <mutex>
#include "../util/log_stream.h"

typedef uint64_t stat_type;

struct Statistics
{

	enum value {
		SEED_HITS, TENTATIVE_MATCHES0, TENTATIVE_MATCHES1, TENTATIVE_MATCHES2, TENTATIVE_MATCHES3, TENTATIVE_MATCHES4, TENTATIVE_MATCHESX, MATCHES, ALIGNED, GAPPED, DUPLICATES,
		GAPPED_HITS, QUERY_SEEDS, QUERY_SEEDS_HIT, REF_SEEDS, REF_SEEDS_HIT, QUERY_SIZE, REF_SIZE, OUT_HITS, OUT_MATCHES, COLLISION_LOOKUPS, QCOV, BIAS_ERRORS, SCORE_TOTAL, ALIGNED_QLEN, PAIRWISE, HIGH_SIM,
		SEARCH_TEMP_SPACE, SECONDARY_HITS, ERASED_HITS, SQUARED_ERROR, CELLS, TARGET_HITS0, TARGET_HITS1, TARGET_HITS2, TARGET_HITS3, TARGET_HITS3_CBS, TARGET_HITS4, TARGET_HITS5, TARGET_HITS6, TIME_GREEDY_EXT, LOW_COMPLEXITY_SEEDS,
		SWIPE_REALIGN, EXT8, EXT16, EXT32, GAPPED_FILTER_TARGETS, GAPPED_FILTER_HITS1, GAPPED_FILTER_HITS2, GROSS_DP_CELLS, NET_DP_CELLS, TIME_TARGET_SORT, TIME_SW, TIME_EXT, TIME_GAPPED_FILTER,
		TIME_LOAD_HIT_TARGETS, TIME_CHAINING, TIME_LOAD_SEED_HITS, TIME_SORT_SEED_HITS, TIME_SORT_TARGETS_BY_SCORE, TIME_TARGET_PARALLEL, TIME_TRACEBACK_SW, TIME_TRACEBACK, HARD_QUERIES, TIME_MATRIX_ADJUST,
		MATRIX_ADJUST_COUNT, MASKED_LAZY, SWIPE_TASKS_TOTAL, SWIPE_TASKS_ASYNC, TRIVIAL_ALN, TIME_EXT_32, EXT_OVERFLOW_8, EXT_WASTED_16, DP_CELLS_8, DP_CELLS_16, DP_CELLS_32, TIME_PROFILE, TIME_ANCHORED_SWIPE,
		TIME_ANCHORED_SWIPE_ALLOC, TIME_ANCHORED_SWIPE_SORT, TIME_ANCHORED_SWIPE_ADD, TIME_ANCHORED_SWIPE_OUTPUT, COUNT
	};

	Statistics()
	{
		reset();
	}

	void reset() {
		std::fill(data_, data_ + COUNT, (stat_type)0);
	}

	Statistics& operator+=(const Statistics &rhs)
	{
		mtx_.lock();
		for(unsigned i=0;i<COUNT;++i)
			data_[i] += rhs.data_[i];
		mtx_.unlock();
		return *this;
	}

	void inc(const value v, stat_type n = 1lu)
	{ data_[v] += n; }

	void max(const value v, stat_type n)
	{
		data_[v] = std::max(data_[v], n);
	}

	stat_type get(const value v) const
	{ return data_[v]; }

	void print() const;

	stat_type data_[COUNT];
	std::mutex mtx_;

};

extern Statistics statistics;
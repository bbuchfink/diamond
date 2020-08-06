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
		SEARCH_TEMP_SPACE, SECONDARY_HITS, ERASED_HITS, SQUARED_ERROR, CELLS, OUTRANKED_HITS, TARGET_HITS0, TARGET_HITS1, TARGET_HITS2, TARGET_HITS3, TARGET_HITS4, TARGET_HITS5, TARGET_HITS6, TIME_GREEDY_EXT, LOW_COMPLEXITY_SEEDS,
		SWIPE_REALIGN, EXT8, EXT16, EXT32, GAPPED_FILTER_TARGETS, GAPPED_FILTER_HITS1, GAPPED_FILTER_HITS2, GROSS_DP_CELLS, NET_DP_CELLS, TIME_TARGET_SORT, TIME_SW, TIME_EXT, TIME_GAPPED_FILTER,
		TIME_LOAD_HIT_TARGETS, TIME_CHAINING, TIME_LOAD_SEED_HITS, TIME_SORT_SEED_HITS, TIME_SORT_TARGETS_BY_SCORE, TIME_TARGET_PARALLEL, TIME_TRACEBACK_SW, TIME_TRACEBACK, COUNT
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

	void print() const
	{
		using std::endl;
		//log_stream << "Used ref size = " << data_[REF_SIZE] << endl;
		//log_stream << "Traceback errors = " << data_[BIAS_ERRORS] << endl;
		//log_stream << "Low complexity seeds  = " << data_[LOW_COMPLEXITY_SEEDS] << endl;
		log_stream << "Hits (filter stage 0) = " << data_[SEED_HITS] << endl;
		log_stream << "Hits (filter stage 1) = " << data_[TENTATIVE_MATCHES1] << " (" << data_[TENTATIVE_MATCHES1]*100.0/ data_[SEED_HITS] << " %)" << endl;
		log_stream << "Hits (filter stage 2) = " << data_[TENTATIVE_MATCHES2] << " (" << data_[TENTATIVE_MATCHES2] * 100.0 / data_[TENTATIVE_MATCHES1] << " %)" << endl;
		log_stream << "Hits (filter stage 3) = " << data_[TENTATIVE_MATCHES3] << " (" << data_[TENTATIVE_MATCHES3] * 100.0 / data_[TENTATIVE_MATCHES2] << " %)" << endl;
		//log_stream << "Hits (filter stage 4) = " << data_[TENTATIVE_MATCHES4] << " (" << data_[TENTATIVE_MATCHES4] * 100.0 / data_[TENTATIVE_MATCHES3] << " %)" << endl;
		log_stream << "Target hits (stage 0) = " << data_[TARGET_HITS0] << endl;
		log_stream << "Target hits (stage 1) = " << data_[TARGET_HITS1] << endl;
		log_stream << "Target hits (stage 2) = " << data_[TARGET_HITS2] << endl;
		log_stream << "Target hits (stage 3) = " << data_[TARGET_HITS3] << endl;
		log_stream << "Target hits (stage 4) = " << data_[TARGET_HITS4] << endl;
		log_stream << "Target hits (stage 5) = " << data_[TARGET_HITS5] << endl;
		log_stream << "Target hits (stage 6) = " << data_[TARGET_HITS6] << endl;
		log_stream << "Swipe realignments    = " << data_[SWIPE_REALIGN] << endl;
		log_stream << "Extensions (8 bit)    = " << data_[EXT8] << endl;
		log_stream << "Extensions (16 bit)   = " << data_[EXT16] << endl;
		log_stream << "Extensions (32 bit)   = " << data_[EXT32] << endl;
#ifdef DP_STAT
		log_stream << "Gross DP Cells        = " << data_[GROSS_DP_CELLS] << endl;
		log_stream << "Net DP Cells          = " << data_[NET_DP_CELLS] << " (" << data_[NET_DP_CELLS] * 100.0 / data_[GROSS_DP_CELLS] << " %)" << endl;
#endif
		log_stream << "Gapped filter (targets) = " << data_[GAPPED_FILTER_TARGETS] << endl;
		log_stream << "Gapped filter (hits) stage 1 = " << data_[GAPPED_FILTER_HITS1] << endl;
		log_stream << "Gapped filter (hits) stage 2 = " << data_[GAPPED_FILTER_HITS2] << endl;
		log_stream << "Time (Load seed hit targets) = " << (double)data_[TIME_LOAD_HIT_TARGETS] / 1e6 << "s (CPU)" << endl;
		log_stream << "Time (Sort targets by score) = " << (double)data_[TIME_SORT_TARGETS_BY_SCORE] / 1e6 << "s (CPU)" << endl;
		log_stream << "Time (Gapped filter)         = " << (double)data_[TIME_GAPPED_FILTER] / 1e6 << "s (CPU)" << endl;
		log_stream << "Time (Chaining)              = " << (double)data_[TIME_CHAINING] / 1e6 << "s (CPU)" << endl;
		log_stream << "Time (DP target sorting)     = " << (double)data_[TIME_TARGET_SORT] / 1e6 << "s (CPU)" << endl;
		log_stream << "Time (Smith Waterman)        = " << (double)data_[TIME_SW] / 1e6 << "s (CPU)" << endl;
		log_stream << "Time (Smith Waterman TB)     = " << (double)data_[TIME_TRACEBACK_SW] / 1e6 << "s (CPU)" << endl;
		log_stream << "Time (Traceback)             = " << (double)data_[TIME_TRACEBACK] / 1e6 << "s (CPU)" << endl;
		log_stream << "Time (Target parallel)       = " << (double)data_[TIME_TARGET_PARALLEL] / 1e6 << "s (wall)" << endl;
		log_stream << "Time (Load seed hits)        = " << (double)data_[TIME_LOAD_SEED_HITS] / 1e6 << "s (wall)" << endl;
		log_stream << "Time (Sort seed hits)        = " << (double)data_[TIME_SORT_SEED_HITS] / 1e6 << "s (wall)" << endl;
		log_stream << "Time (Extension)             = " << (double)data_[TIME_EXT] / 1e6 << "s (wall)" << endl;
		//log_stream << "Time (greedy extension)      = " << data_[TIME_GREEDY_EXT]/1e9 << "s" << endl;
		//log_stream << "Gapped hits = " << data_[GAPPED_HITS] << endl;
		//log_stream << "Overlap hits = " << data_[DUPLICATES] << endl;
		//log_stream << "Secondary hits = " << data_[SECONDARY_HITS] << endl;
		//log_stream << "Erased hits = " << data_[ERASED_HITS] << endl;
		//log_stream << "High similarity hits = " << data_[HIGH_SIM] << endl;
		//log_stream << "Net hits = " << data_[OUT_HITS] << endl;
		//log_stream << "Matches = " << data_[OUT_MATCHES] << endl;
		//log_stream << "Total score = " << data_[SCORE_TOTAL] << endl;
		//log_stream << "Aligned query len = " << data_[ALIGNED_QLEN] << endl;
		//log_stream << "Gapped matches = " << data_[GAPPED] << endl;
		//log_stream << "MSE = " << (double)data_[SQUARED_ERROR] / (double)data_[OUT_HITS] << endl;
		//log_stream << "Cells = " << data_[CELLS] << endl;
		verbose_stream << "Temporary disk space used (search): " << (double)data_[SEARCH_TEMP_SPACE] / (1 << 30) << " GB" << endl;
		if(data_[OUTRANKED_HITS])
			log_stream << "Outranked hits = " << data_[OUTRANKED_HITS] << " (" << data_[OUTRANKED_HITS]*100.0/ data_[PAIRWISE] << "%)" << endl;
		message_stream << "Reported " << data_[PAIRWISE] << " pairwise alignments, " << data_[MATCHES] << " HSPs." << endl;
		message_stream << data_[ALIGNED] << " queries aligned." << endl;
	}

	stat_type data_[COUNT];
	std::mutex mtx_;

};

extern Statistics statistics;
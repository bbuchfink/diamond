/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#ifndef STATISTICS_H_
#define STATISTICS_H_

#include <algorithm>
#include <stdint.h>
#include <string.h>
#include "../util/tinythread.h"
#include "../util/log_stream.h"

typedef uint64_t stat_type;

struct Statistics
{

	enum value {
		SEED_HITS, TENTATIVE_MATCHES0, TENTATIVE_MATCHES1, TENTATIVE_MATCHES2, TENTATIVE_MATCHES3, TENTATIVE_MATCHES4, TENTATIVE_MATCHESX, MATCHES, ALIGNED, GAPPED, DUPLICATES,
		GAPPED_HITS, QUERY_SEEDS, QUERY_SEEDS_HIT, REF_SEEDS, REF_SEEDS_HIT, QUERY_SIZE, REF_SIZE, OUT_HITS, OUT_MATCHES, COLLISION_LOOKUPS, QCOV, BIAS_ERRORS, SCORE_TOTAL, ALIGNED_QLEN, PAIRWISE, HIGH_SIM,
		TEMP_SPACE, SECONDARY_HITS, ERASED_HITS, SQUARED_ERROR, CELLS, OUTRANKED_HITS, TARGET_HITS0, TARGET_HITS1, TARGET_HITS2, TIME_GREEDY_EXT, COUNT
	};

	Statistics()
	{ memset(data_, 0, sizeof(data_)); }

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
		//log_stream << "Used ref size = " << data_[REF_SIZE] << endl;
		//log_stream << "Traceback errors = " << data_[BIAS_ERRORS] << endl;
		log_stream << "Hits (filter stage 0) = " << data_[SEED_HITS] << endl;
		log_stream << "Hits (filter stage 1) = " << data_[TENTATIVE_MATCHES1] << " (" << data_[TENTATIVE_MATCHES1]*100.0/ data_[SEED_HITS] << " %)" << endl;
		log_stream << "Hits (filter stage 2) = " << data_[TENTATIVE_MATCHES2] << " (" << data_[TENTATIVE_MATCHES2] * 100.0 / data_[TENTATIVE_MATCHES1] << " %)" << endl;
		log_stream << "Hits (filter stage x) = " << data_[TENTATIVE_MATCHESX] << " (" << data_[TENTATIVE_MATCHESX] * 100.0 / data_[TENTATIVE_MATCHES2] << " %)" << endl;
		log_stream << "Hits (filter stage 3) = " << data_[TENTATIVE_MATCHES3] << " (" << data_[TENTATIVE_MATCHES3] * 100.0 / data_[TENTATIVE_MATCHESX] << " %)" << endl;
		log_stream << "Hits (filter stage 4) = " << data_[TENTATIVE_MATCHES4] << " (" << data_[TENTATIVE_MATCHES4] * 100.0 / data_[TENTATIVE_MATCHES3] << " %)" << endl;
		log_stream << "Target hits (stage 0) = " << data_[TARGET_HITS0] << endl;
		log_stream << "Target hits (stage 1) = " << data_[TARGET_HITS1] << endl;
		log_stream << "Target hits (stage 2) = " << data_[TARGET_HITS2] << endl;
		log_stream << "Time (greedy extension) = " << data_[TIME_GREEDY_EXT]/1e9 << "s" << endl;
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
		log_stream << "MSE = " << (double)data_[SQUARED_ERROR] / (double)data_[OUT_HITS] << endl;
		//log_stream << "Cells = " << data_[CELLS] << endl;
		verbose_stream << "Temporary disk space used: " << (double)data_[TEMP_SPACE] / (1 << 30) << " GB" << endl;
		log_stream << "Outranked hits = " << data_[OUTRANKED_HITS] << " (" << data_[OUTRANKED_HITS]*100.0/ data_[PAIRWISE] << "%)" << endl;
		message_stream << "Reported " << data_[PAIRWISE] << " pairwise alignments, " << data_[MATCHES] << " HSPs." << endl;
		message_stream << data_[ALIGNED] << " queries aligned." << endl;
	}

	stat_type data_[COUNT];
	tthread::mutex mtx_;

};

extern Statistics statistics;

#endif /* STATISTICS_H_ */

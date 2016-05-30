/****
Copyright (c) 2014, University of Tuebingen
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****
Author: Benjamin Buchfink
****/

#ifndef STATISTICS_H_
#define STATISTICS_H_

#include "../util/tinythread.h"

typedef uint64_t stat_type;

struct Statistics
{

	enum value { SEED_HITS, TENTATIVE_MATCHES0, TENTATIVE_MATCHES1, TENTATIVE_MATCHES2, TENTATIVE_MATCHES3, TENTATIVE_MATCHES4, MATCHES, ALIGNED, GAPPED, DUPLICATES,
		GAPPED_HITS, QUERY_SEEDS, QUERY_SEEDS_HIT, REF_SEEDS, REF_SEEDS_HIT, QUERY_SIZE, REF_SIZE, OUT_HITS, OUT_MATCHES, COLLISION_LOOKUPS, QCOV, BIAS_ERRORS, SCORE_TOTAL, ALIGNED_QLEN, PAIRWISE, COUNT };

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

	stat_type get(const value v) const
	{ return data_[v]; }

	void print() const
	{
		log_stream << "Used ref size = " << data_[REF_SIZE] << endl;
		log_stream << "Traceback errors = " << data_[BIAS_ERRORS] << endl;
		verbose_stream << "Hits (filter stage 0) = " << data_[SEED_HITS] << endl;
		verbose_stream << "Hits (filter stage 1) = " << data_[TENTATIVE_MATCHES1] << " (" << data_[TENTATIVE_MATCHES1]*100.0/ data_[SEED_HITS] << " %)" << endl;
		verbose_stream << "Hits (filter stage 2) = " << data_[TENTATIVE_MATCHES2] << " (" << data_[TENTATIVE_MATCHES2] * 100.0 / data_[TENTATIVE_MATCHES1] << " %)" << endl;
		verbose_stream << "Hits (filter stage 3) = " << data_[TENTATIVE_MATCHES3] << " (" << data_[TENTATIVE_MATCHES3] * 100.0 / data_[TENTATIVE_MATCHES2] << " %)" << endl;
		verbose_stream << "Hits (filter stage 4) = " << data_[TENTATIVE_MATCHES4] << " (" << data_[TENTATIVE_MATCHES4] * 100.0 / data_[TENTATIVE_MATCHES3] << " %)" << endl;
		log_stream << "Gapped hits = " << data_[GAPPED_HITS] << endl;
		log_stream << "Overlap hits = " << data_[DUPLICATES] << endl;
		log_stream << "Net hits = " << data_[OUT_HITS] << endl;
		log_stream << "Matches = " << data_[OUT_MATCHES] << endl;
		log_stream << "Total score = " << data_[SCORE_TOTAL] << endl;
		log_stream << "Aligned query len = " << data_[ALIGNED_QLEN] << endl;
		log_stream << "Gapped matches = " << data_[GAPPED] << endl;
		message_stream << "Reported " << data_[PAIRWISE] << " pairwise alignments, " << data_[MATCHES] << " HSSPs." << endl;
		message_stream << data_[ALIGNED] << " queries aligned." << endl;
	}

	stat_type data_[COUNT];
	tthread::mutex mtx_;

};

extern Statistics statistics;

#endif /* STATISTICS_H_ */

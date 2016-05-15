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

#ifndef COMPLEXITY_FILTER_H_
#define COMPLEXITY_FILTER_H_

#include "../blast/blast_seg.h"
#include "../blast/blast_filter.h"
#include "../basic/value.h"
#include "../data/sequence_set.h"
#include "thread.h"

struct Complexity_filter
{

	Complexity_filter()
	{ blast_seg_ = SegParametersNewAa(); }

	~Complexity_filter()
	{ SegParametersFree(blast_seg_); }

	unsigned filter(Letter *seq, uint32_t length) const
	{
		BlastSeqLoc *seg_locs;
		SeqBufferSeg ((uint8_t*) seq, length, 0u, blast_seg_, &seg_locs);
		unsigned nMasked = 0;

		if(seg_locs) {
			BlastSeqLoc *l = seg_locs;
			do {
				for(signed i=l->ssr->left;i<=l->ssr->right;i++) {
					nMasked++;
					seq[i] = value_traits.mask_char;
				}
			} while((l=l->next) != 0);
			BlastSeqLocFree(seg_locs);
		}
		return nMasked;
	}

	static const Complexity_filter& get()
	{ return instance; }

	void run(Sequence_set &seqs) const
	{
		Filter_context context (seqs, *this);
		launch_scheduled_thread_pool(context, (unsigned)seqs.get_length(), config.threads_);
	}

private:

	struct Filter_context
	{
		Filter_context(Sequence_set &seqs, const Complexity_filter &filter):
			seqs (seqs),
			filter (filter)
		{ }
		void operator()(unsigned thread_id, unsigned i)
		{
			filter.filter(seqs.ptr(i), (unsigned)seqs.length(i));
		}
		Sequence_set &seqs;
		const Complexity_filter &filter;
	};

	SegParameters *blast_seg_;

	static const Complexity_filter instance;

};

#endif /* COMPLEXITY_FILTER_H_ */

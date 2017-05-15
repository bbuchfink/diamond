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

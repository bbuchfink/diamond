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

#include "../algo/blast/core/blast_seg.h"
#include "../algo/blast/core/blast_filter.h"
#include "../basic/value.h"

template<class _val>
struct Complexity_filter
{
	unsigned filter(vector<_val> &seq) const
	{ return 0; }
	static const Complexity_filter& get()
	{ return instance; }
	void run(String_set<_val> &seqs) const
	{ }
	static const Complexity_filter instance;
};

template<>
struct Complexity_filter<Amino_acid>
{

	Complexity_filter()
	{ blast_seg_ = SegParametersNewAa(); }

	~Complexity_filter()
	{ SegParametersFree(blast_seg_); }

	unsigned filter(sequence<Amino_acid> seq) const
	{
		BlastSeqLoc *seg_locs;
		SeqBufferSeg ((uint8_t*) seq.data(), seq.length(), 0, blast_seg_, &seg_locs);
		unsigned nMasked = 0;

		if(seg_locs) {
			BlastSeqLoc *l = seg_locs;
			do {
				for(signed i=l->ssr->left;i<=l->ssr->right;i++) {
					nMasked++;
					seq[i] = Value_traits<Amino_acid>::MASK_CHAR;
				}
			} while((l=l->next) != 0);
			BlastSeqLocFree(seg_locs);
		}
		return nMasked;
	}

	static const Complexity_filter& get()
	{ return instance; }

	void run(String_set<Amino_acid> &seqs) const
	{
		Filter_context context (seqs, *this);
		launch_scheduled_thread_pool(context, seqs.get_length(), program_options::threads());
	}

private:

	struct Filter_context
	{
		Filter_context(String_set<Amino_acid> &seqs, const Complexity_filter &filter):
			seqs (seqs),
			filter (filter)
		{ }
		void operator()(unsigned thread_id, unsigned i)
		{
			filter.filter(seqs[i]);
		}
		String_set<Amino_acid> &seqs;
		const Complexity_filter &filter;
	};

	SegParameters *blast_seg_;

	static const Complexity_filter instance;

};

const Complexity_filter<Amino_acid> Complexity_filter<Amino_acid>::instance;
#ifdef NDEBUG
template<> const Complexity_filter<Nucleotide> Complexity_filter<Nucleotide>::instance;
#else
template<typename _val> const Complexity_filter<_val> Complexity_filter<_val>::instance;
#endif


#endif /* COMPLEXITY_FILTER_H_ */

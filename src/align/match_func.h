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

#ifndef MATCH_FUNC_H_
#define MATCH_FUNC_H_

#include "../basic/match.h"

int blast_frame(unsigned frame)
{ return frame <= 2 ? (int)frame+1 : 2-(int)frame; }

template<typename _val>
void anchored_transform(local_match<_val> &l, unsigned subject_pos, unsigned query_pos)
{
	l.query_begin_ = query_pos - l.query_begin_;
	l.subject_begin_ = subject_pos - l.subject_begin_;
}

template<typename _val>
void to_source_space(local_match<_val> &l, unsigned frame, unsigned source_len)
{
	if(frame == 1) {
		l.query_begin_ = source_len - l.query_begin_ - 1;
		l.query_len_ *= -1;
	}
}

template<>
void to_source_space<Amino_acid>(local_match<Amino_acid> &l, unsigned frame, unsigned dna_len)
{
	if(!query_translated())
		return;
	int query_begin_dna; //, query_end_dna;
	signed f = frame <= 2 ? frame+1 : 2-frame;
	if (f > 0) {
		query_begin_dna = (f-1) + 3 * l.query_begin_;
		//query_end_dna = (f-1) + 3 * (l.query_begin_+l.query_len_-1) + 3;
	} else {
		query_begin_dna = dna_len + f - 3 * l.query_begin_;
		//query_end_dna = dna_len + (f + 1) - 3 * (l.query_begin_+l.query_len_-1) - 2;
		l.query_len_ *= -1;
	}
	l.query_begin_ = query_begin_dna;
	l.query_len_*= 3;
}

template<typename _val>
unsigned query_translated_begin(unsigned query_begin, unsigned frame, unsigned dna_len)
{
	if(frame == 0)
		return query_begin;
	else
		return dna_len-query_begin-1;
}

template<>
unsigned query_translated_begin<Amino_acid>(unsigned query_begin, unsigned frame, unsigned dna_len)
{
	if(!query_translated())
		return query_begin;
	int f = frame <= 2 ? frame+1 : 2-frame;
	if (f > 0)
		return (query_begin - (f-1))/3;
	else
		return (dna_len + f - query_begin)/3;
}

#endif /* MATCH_FUNC_H_ */

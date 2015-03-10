/****
Copyright (c) 2014, University of Tuebingen
Author: Benjamin Buchfink
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/

#ifndef SCORE_TRAITS_H_
#define SCORE_TRAITS_H_

template<typename _val>
uint8_t blast_seq_code()
{ return 0; }

template<>
uint8_t blast_seq_code<Amino_acid>()
{ return BLASTAA_SEQ_CODE; }

template<typename _val>
int16_t blast_load_karlin_blk(Blast_KarlinBlk* kbp,
		Blast_KarlinBlk* kbp_ungap,
		int gap_open,
		int gap_extend,
		int reward,
		int penalty,
		const char *matrix)
{ return 0; }

template<>
int16_t blast_load_karlin_blk<Amino_acid>(Blast_KarlinBlk* kbp,
		Blast_KarlinBlk* kbp_ungap,
		int gap_open,
		int gap_extend,
		int reward,
		int penalty,
		const char *matrix)
{
	return Blast_KarlinBlkGappedLoadFromTables(kbp,
			gap_open,
			gap_extend,
			matrix);
}

template<typename _val>
const uint8_t* blast_alphabet()
{ return 0; }

template<>
const uint8_t* blast_alphabet<Amino_acid>()
{ return AMINOACID_TO_NCBISTDAA; }

#ifdef EXTRA
#include "../../extra/score_traits.h"
#endif

#endif /* SCORE_TRAITS_H_ */

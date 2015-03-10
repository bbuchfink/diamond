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

#ifndef ALIGN_H_
#define ALIGN_H_

#include "../data/reference.h"
#include "../basic/statistics.h"
#include "../basic/score_matrix.h"
#include "../basic/shape_config.h"
#include "../search/sse_dist.h"
#include "../search/collision.h"
#include "../search/hit_filter.h"
#include "../search/align_ungapped.h"

template<typename _val, typename _locr, typename _locq, typename _locl>
void align(const _locq q_pos,
	  const _val *query,
	  _locr s,
	  Statistics &stats,
	  const unsigned sid,
	  hit_filter<_val,_locr,_locq,_locl> &hf)
{
	stats.inc(Statistics::TENTATIVE_MATCHES0);
	const _val* subject = ref_seqs<_val>::data_->data(s);

	if(fast_match(query, subject) < program_options::min_identities)
		return;

	stats.inc(Statistics::TENTATIVE_MATCHES1);

	unsigned delta, len;
	int score;
	if((score = xdrop_ungapped<_val,_locr,_locq>(query, subject, shape_config::get().get_shape(sid).length_, delta, len)) < program_options::min_ungapped_raw_score)
		return;

	if(!is_primary_hit<_val,_locr>(query-delta, subject-delta, delta, sid, len))
		return;

	stats.inc(Statistics::TENTATIVE_MATCHES2);
	hf.push(s, score);
}

#endif /* ALIGN_H_ */

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

#ifndef FILTER_HIT_H_
#define FILTER_HIT_H_

#include "../data/reference.h"
#include "../basic/statistics.h"
#include "../basic/score_matrix.h"
#include "../basic/shape_config.h"
#include "../search/sse_dist.h"
#include "../search/collision.h"
#include "../search/hit_filter.h"
#include "../dp/dp.h"

inline void align(Loc q_pos,
	  const Letter *query,
	  Loc s,
	  Statistics &stats,
	  const unsigned sid,
	  hit_filter &hf)
{
	const Letter* subject = ref_seqs::data_->data(s);

	if(fast_match(query, subject) < config.min_identities)
		return;

	stats.inc(Statistics::TENTATIVE_MATCHES1);

	unsigned delta, len;
	int score;
	if((score = xdrop_ungapped(query, subject, shapes[sid].length_, delta, len)) < config.min_ungapped_raw_score)
		return;

	stats.inc(Statistics::TENTATIVE_MATCHES2);

	if(!is_primary_hit(query-delta, subject-delta, delta, sid, len))
		return;

	stats.inc(Statistics::TENTATIVE_MATCHES3);
	hf.push(s, score);
}

#endif

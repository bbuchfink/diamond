/****
Copyright (c) 2014-2015, University of Tuebingen
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

#ifndef LINK_SEGMENTS_H_
#define LINK_SEGMENTS_H_

#include <vector>
#include "../util/map.h"
#include "../basic/match.h"

using std::vector;

static const size_t MAX_LINKING_OVERLAP = 10;

/**double link_segments(Segment &h1, Segment &h2)
{
	if(h1.strand() == h2.strand() && h2.top_evalue_ == -1) {
		Segment *p = &h1;
		while(p != 0) {
			if(p->query_range().overlap(h2.query_range())/3 >= 1
				|| p->subject_range().overlap(h2.subject_range()) > MAX_LINKING_OVERLAP)
				return std::numeric_limits<double>::max();
			p = p->next_;
		}
		double ev = h1.evalue_ * h2.evalue_;
		h2.next_ = h1.next_;
		h1.next_ = &h2;
		p = &h1;
		while(p != 0) {
			p->evalue_ = ev;
			p->top_evalue_ = ev;
			p = p->next_;
		}
		return ev;
	}
	return std::numeric_limits<double>::max();
}*/

inline void link_segments(const vector<Segment>::iterator &begin, const vector<Segment>::iterator &end)
{
	int max_score = begin->score_;
	/*for(typename vector<match<_val> >::iterator i=begin; i<end; ++i)
		if(i->top_evalue_ == -1)
			for(typename vector<match<_val> >::iterator j=i+1; j<end; ++j)
				min_ev = std::min(min_ev, link_segments(*i, *j));*/
	for(vector<Segment>::iterator i=begin; i<end; ++i)
		i->top_score_ = max_score;
}

inline void link_segments(vector<Segment> &hsp_list)
{
	typedef Map<vector<Segment>::iterator,Segment::Subject> Hsp_map;
	std::sort(hsp_list.begin(), hsp_list.end(), Segment::comp_subject);
	Hsp_map hsp_map (hsp_list.begin(), hsp_list.end());
	Hsp_map::Iterator it = hsp_map.begin();
	while(it.valid()) {
		link_segments(it.begin(), it.end());
		++it;
	}
}

#endif /* LINK_SEGMENTS_H_ */

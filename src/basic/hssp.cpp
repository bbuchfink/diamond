/****
Copyright (c) 2016, Benjamin Buchfink
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

#include "match.h"

bool local_match::pass_through(const Diagonal_segment &d, const vector<char> &transcript_buf)
{
	if (intersect(d.query_range(), query_range(FORWARD)).length() != d.len
		|| intersect(d.subject_range(), subject_range()).length() != d.len)
		return false;

	Link_iterator it(transcript_right_, transcript_left_, transcript_buf);
	const int subject_end = d.subject_pos + d.len;
	const int diag = d.diag();
	int i = query_begin_, j = subject_begin_;
	while (true) {
		if (j >= d.subject_pos) {
			if (j >= subject_end)
				return true;
			if (j-i != diag)
				return false;
		}
		++it;
		if (!it.good())
			break;
		switch (*it) {
		case op_match:
			++i;
			++j;
			break;
		case op_deletion:
			++j;
			break;
		case op_insertion:
			++i;
		}		
	}
	return true;
}

bool local_match::is_weakly_enveloped(const local_match &j)
{
	static const double overlap_factor = 0.9;
	return score_ <= j.score_
		&& subject_range().overlap_factor(j.subject_range()) >= overlap_factor
		&& query_range(FORWARD).overlap_factor(j.query_range(FORWARD)) >= overlap_factor;
}

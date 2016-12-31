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

#ifndef ALIGN_STRUCT_H_
#define ALIGN_STRUCT_H_

#include "../basic/match.h"

struct local_match : public Hsp_data
{
	local_match() :
		query_anchor_(0),
		subject_(0)
	{ }
	local_match(int score) :
		Hsp_data(score)
	{ }
	local_match(int query_anchor, int subject_anchor, const Letter *subject, unsigned total_subject_len = 0) :
		total_subject_len_(total_subject_len),
		query_anchor_(query_anchor),
		subject_anchor(subject_anchor),
		subject_(subject)
	{ }
	local_match(unsigned len, unsigned query_begin, unsigned query_len, unsigned subject_len, unsigned gap_openings, unsigned identities, unsigned mismatches, signed subject_begin, signed score) :
		Hsp_data(score),
		query_anchor_(0),
		subject_(0)
	{ }
	unsigned total_subject_len_;
	signed query_anchor_, subject_anchor;
	const Letter *subject_;
};

#endif
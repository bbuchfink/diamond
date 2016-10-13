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

#include "dp.h"
#include "../basic/score_matrix.h"

int xdrop_ungapped(const Letter *query, const Letter *subject, unsigned seed_len, unsigned &delta, unsigned &len)
{
	int score(0), st(0);
	unsigned n(0);
	delta = 0;

	const Letter *q(query - 1), *s(subject - 1);
	const unsigned window_left = std::max(config.window, (unsigned)config.seed_anchor) - config.seed_anchor;
	while (score - st < config.raw_ungapped_xdrop
		&& delta < window_left
		&& *q != '\xff'
		&& *s != '\xff')
	{
		st += score_matrix(*q, *s);
		score = std::max(score, st);
		--q;
		--s;
		++delta;
	}

	q = query + seed_len;
	s = subject + seed_len;
	st = score;
	assert(seed_len >= config.seed_anchor);
	const unsigned window_right = std::max(config.window, seed_len - config.seed_anchor) - (seed_len - config.seed_anchor);
	while (score - st < config.raw_ungapped_xdrop
		&& n < window_right
		&& *q != '\xff'
		&& *s != '\xff')
	{
		st += score_matrix(*q, *s);
		score = std::max(score, st);
		++q;
		++s;
		++n;
	}

	for (unsigned i = 0; i<seed_len; ++i)
		score += score_matrix(query[i], subject[i]);

	len = delta + n + seed_len;
	return score;
}

int xdrop_ungapped(const Letter *query, const Letter *subject, unsigned &delta, unsigned &len)
{
	int score(0), st(0), n=1;
	delta = 0;
	
	const Letter *q(query - 1), *s(subject - 1);
	while (score - st < config.raw_ungapped_xdrop
		&& *q != '\xff'
		&& *s != '\xff')
	{
		st += score_matrix(*q, *s);
		if (st > score) {
			score = st;
			delta = n;
		}
		--q;
		--s;
		++n;
	}

	q = query;
	s = subject;
	st = score;
	n = 1;
	len = 0;
	while (score - st < config.raw_ungapped_xdrop
		&& *q != '\xff'
		&& *s != '\xff')
	{
		st += score_matrix(*q, *s);
		if (st > score) {
			score = st;
			len = n;
		}
		++q;
		++s;
		++n;
	}
	len += delta;
	return score;
}

int xdrop_ungapped_right(const Letter *query, const Letter *subject, int &len)
{
	int score(0), st(0), n = 1;
	len = 0;
	
	const Letter *q = query;
	const Letter *s = subject;
	
	while (score - st < config.raw_ungapped_xdrop
		&& *q != '\xff'
		&& *s != '\xff')
	{
		st += score_matrix(*q, *s);
		if (st > score) {
			score = st;
			len = n;
		}
		++q;
		++s;
		++n;
	}
	return score;
}
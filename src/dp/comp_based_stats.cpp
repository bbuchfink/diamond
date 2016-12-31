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

double background_scores[20];
const double background_freq[] = { 0.0844581,0.0581912,0.0421072,0.0546748,0.0146359,0.040118,0.0621211,0.0669379,0.0225159,0.0547866,0.0957934,0.0523275,0.0218629,0.038769,0.0505311,
0.0760908,0.0573267,0.0127314,0.0295317,0.0644889 };

void init_cbs()
{	
	for (unsigned i = 0; i < 20; ++i) {
		background_scores[i] = 0;
		for (unsigned j = 0; j < 20; ++j)
			background_scores[i] += background_freq[j] * score_matrix(i, j);
	}
}

struct Vector_scores
{
	Vector_scores()
	{
		memset(scores, 0, sizeof(scores));
	}
	Vector_scores& operator+=(Letter l)
	{
		for (unsigned i = 0; i < 20; ++i)
			scores[i] += score_matrix(l, i);
		return *this;
	}
	Vector_scores& operator-=(Letter l)
	{
		for (unsigned i = 0; i < 20; ++i)
			scores[i] -= score_matrix(l, i);
		return *this;
	}
	int scores[20];
};

Bias_correction::Bias_correction(const sequence &seq):
	vector<float>(seq.length())
{
	Vector_scores scores;
	const unsigned window_half = std::min(20u, (unsigned)seq.length() - 1);
	unsigned n = 0;
	unsigned h = 0, m = 0, t = 0, l = (unsigned)seq.length();
	while (n < window_half && h < l) {
		++n;
		scores += seq[h];
		++h;
	}
	while (n < 41 && h < l) {
		++n;
		scores += seq[h];
		const Letter r = seq[m];
		if (r < 20)
			this->operator[](m) = (float)background_scores[(int)r] - float(scores.scores[(int)r] - score_matrix(r, r)) / (n - 1);
		++h;
		++m;
	}
	while (h < l) {
		scores += seq[h];
		scores -= seq[t];
		const Letter r = seq[m];
		if (r < 20)
			this->operator[](m) = (float)background_scores[(int)r] - float(scores.scores[(int)r] - score_matrix(r, r)) / (n - 1);
		++h;
		++t;
		++m;
	}
	while (m < l && n>21) {
		--n;
		scores -= seq[t];
		const Letter r = seq[m];
		if (r < 20)
			this->operator[](m) = (float)background_scores[(int)r] - float(scores.scores[(int)r] - score_matrix(r, r)) / (n - 1);
		++t;
		++m;
	}
	while (m < l) {
		const Letter r = seq[m];
		if (r < 20)
			this->operator[](m) = (float)background_scores[(int)r] - float(scores.scores[(int)r] - score_matrix(r, r)) / (n - 1);
		++m;
	}
}

int Bias_correction::operator()(const Hsp_data &hsp) const
{
	float s = 0;
	for (Hsp_data::Iterator i = hsp.begin(); i.good(); ++i) {
		switch (i.op()) {
		case op_match:
		case op_substitution:
			s += (*this)[i.query_pos];
		default:
			;
		}
	}
	return (int)s;
}
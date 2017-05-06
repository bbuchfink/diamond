/****
Copyright (c) 2017, Benjamin Buchfink
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

using std::pair;

struct Banded_dp_matrix
{

	struct Column_iterator
	{
		inline Column_iterator(const pair<int*, int*> &score, int *hgap) :
			score_(score),
			hgap_(hgap + 1)
		{
		}
		inline int& score()
		{
			return *score_.second;
		}
		inline int diag() const
		{
			return *score_.first;
		}
		inline int hgap_in() const
		{
			return *hgap_;
		}
		inline int& hgap_out()
		{
			return *(hgap_ - 1);
		}
		inline void operator++()
		{
			++score_.first;
			++score_.second;
			++hgap_;
		}
	private:
		pair<int*, int*> score_;
		int *hgap_;
	};

	inline Column_iterator column(int j, int offset)
	{
		return Column_iterator(std::make_pair(&score_[j*band_ + offset], &score_[(j + 1)*band_ + offset]), hgap_.data() + offset);
	}

	inline Banded_dp_matrix(int band, int cols) :
		band_(band),
		score_(band*(cols + 1)),
		hgap_(band + 1)
	{}

private:

	const int band_;
	vector<int> score_, hgap_;

};

void banded_sw(const sequence &query, const sequence &subject, int d_begin, int d_end, int j_begin, int j_end, Hsp_data &out)
{
	using std::max;
	assert(d_end > d_begin);
	const int slen = (int)subject.length(),
		qlen = (int)query.length(),
		band = d_end - d_begin,
		gap_open = score_matrix.gap_open() + score_matrix.gap_extend(),
		gap_extend = score_matrix.gap_extend();
	int i0 = d_begin + j_begin, score = 0;
	Banded_dp_matrix mtx(band, slen);
	for (int j = j_begin; j < j_end; ++j, ++i0) {
		const int i1 = std::min(i0 + band, qlen);
		int i = std::max(i0, 0), vgap = 0;
		Banded_dp_matrix::Column_iterator it = mtx.column(j, i - i0);
		for (; i < i1; ++i, ++it) {
			const int match_score = score_matrix(query[i], subject[j]);
			const int hgap = it.hgap_in();
			const int s = max(max(max(it.diag() + match_score, vgap), hgap), 0);
			const int open = s - gap_open;
			vgap = max(vgap - gap_extend, open);
			it.hgap_out() = max(hgap - gap_extend, open);
			it.score() = s;
			score = max(score, s);
		}
	}
	out.score = score;
}
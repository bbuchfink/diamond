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

struct Banded_traceback_matrix
{
	Banded_traceback_matrix(const vector<int> &data, int band, int i0) :
		data_(data),
		band_(band),
		i0_(i0)
	{ }
	int operator()(int i, int j) const
	{
		return data_[(j + 1)*band_ + (i - (i0_ + j))];
	}
	bool in_band(int i, int j) const
	{
		return i >= 0 && j >= 0 && i >= i0_ + j && i < i0_ + j + band_;
	}
	void print(int col, int row) const
	{
		for (unsigned j = 0; j <= row; ++j) {
			for (unsigned i = 0; i <= col; ++i)
				printf("%4i", in_band(i, j) ? this->operator()(i, j) : 0);
			printf("\n");
		}
	}
private:
	const vector<int> &data_;
	const int band_, i0_;
};

int have_gap(const Banded_traceback_matrix &dp,
	int i,
	int j,
	int &l)
{
	const int score = dp(i, j);
	l = 1;
	int j1 = j - 1, i1 = i - 1, g = score_matrix.gap_open() + score_matrix.gap_extend();
	bool bv, bh;
	while ((bv = dp.in_band(i1, j)) || (bh = dp.in_band(i, j1))) {
		if (bv && score == dp(i1, j) - g)
			return 0;
		else if (bh && score == dp(i, j1) - g)
			return 1;
		--j1;
		--i1;
		++l;
		++g;
	}
	return -1;
}

void traceback(const sequence &query,
	const sequence &subject,
	const vector<int> &scores,
	int band,
	int i0,
	int i,
	int j,
	Hsp_data &l)
{
	Banded_traceback_matrix dp(scores, band, i0);
	//dp.print(i, j);
	l.query_range.end_ = i + 1;
	l.subject_range.end_ = j + 1;
	l.transcript.clear();

	int gap_len, score;

	while ((score = dp(i, j)) > 0) {
		const int match_score = score_matrix(query[i], subject[j]);
		//printf("i=%i j=%i score=%i subject=%c query=%c\n",i,j,dp(i, j),Value_traits<_val>::ALPHABET[ls],Value_traits<_val>::ALPHABET[lq]);

		if (score == match_score + dp(i - 1, j - 1)) {
			if (query[i] == subject[j]) {
				l.transcript.push_back(op_match);
				++l.identities;
				++l.positives;
			}
			else {
				l.transcript.push_back(op_substitution, subject[j]);
				++l.mismatches;
				if (match_score > 0)
					++l.positives;
			}
			--i;
			--j;
			++l.length;
		}
		else {
			const int g = have_gap(dp, i, j, gap_len);
			if (g == -1)
				throw std::runtime_error("Traceback error.");
			++l.gap_openings;
			l.length += gap_len;
			l.gaps += gap_len;
			if (g == 0) {
				i -= gap_len;
				l.transcript.push_back(op_insertion, (unsigned)gap_len);
			}
			else {
				for (; gap_len > 0; gap_len--)
					j--;
					l.transcript.push_back(op_deletion, subject[j--]);
			}
		}
	}

	l.query_range.begin_ = i + 1;
	l.subject_range.begin_ = j + 1;
	l.transcript.reverse();
	l.transcript.push_terminator();
}

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

	const vector<int>& scores() const
	{
		return score_;
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
		qlen = (int)query.length();
	d_begin = std::max(d_begin, -(slen - 1));
	d_end = std::min(d_end, qlen);
	const int band = d_end - d_begin,
		gap_open = score_matrix.gap_open() + score_matrix.gap_extend(),
		gap_extend = score_matrix.gap_extend();
	int i0 = d_begin + j_begin, score = 0, max_i, max_j;
	int n=0;
	Banded_dp_matrix mtx(band, slen);
	for (int j = j_begin; j < j_end; ++j, ++i0) {
		const int i1 = std::min(i0 + band, qlen);
		int i = std::max(i0, 0), vgap = 0;
		Banded_dp_matrix::Column_iterator it = mtx.column(j, i - i0);
		for (; i < i1; ++i, ++it) {
			const int match_score = score_matrix(query[i], subject[j]);
			int hgap = it.hgap_in();
			//const int s = max(max(max(it.diag() + match_score, vgap), hgap), 0);
			int s = it.diag() + match_score;
			if (s < hgap)
				s = hgap;
			if (s < vgap)
				s = vgap;			
			if (s < 0)
				s = 0;
			const int open = s - gap_open;
			//vgap = max(vgap - gap_extend, open);
			vgap -= gap_extend;
			if (vgap < open)
				vgap = open;
			//it.hgap_out() = max(hgap - gap_extend, open);
			hgap -= gap_extend;
			if (hgap < open)
				hgap = open;
			it.hgap_out() = hgap;
			it.score() = s;
			//score = max(score, s);
			if (s > score) {
				score = s;
				max_i = i;
				max_j = j;
			}
		}
	}
	out.score = score;
	traceback(query, subject, mtx.scores(), band, d_begin + j_begin, max_i, max_j, out);
}
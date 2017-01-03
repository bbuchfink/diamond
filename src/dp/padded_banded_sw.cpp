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

struct Padded_banded_DP_matrix
{

	struct Column_iterator
	{
		Column_iterator(unsigned column, int* hgap_front, int* score_front, unsigned row_pos, unsigned row_end, unsigned delta) :
			row_pos_(row_pos),
			row_end_(row_end),
			delta_(delta),
			hgap_ptr_(hgap_front),
			score_ptr_(score_front),
			d_(delta > 0 ? *score_front : 0)
		{ }
		inline bool at_end() const
		{
			return row_pos_ >= row_end_;
		}
		inline void operator++()
		{
			++row_pos_; ++hgap_ptr_; ++score_ptr_;
		}
		inline int hgap() const
		{
			return *(hgap_ptr_ + delta_);
		}
		inline int diag() const
		{
			return d_;
		}
		inline void set_hgap(int x)
		{
			*hgap_ptr_ = x;
		}
		inline void set_score(int x)
		{
			d_ = *(score_ptr_ + delta_); *score_ptr_ = x;
		}
		unsigned row_pos_, row_end_, delta_;
		int *hgap_ptr_, *score_ptr_, d_;
	};

	Padded_banded_DP_matrix(unsigned columns, unsigned rows, unsigned band, unsigned padding) :
		rows_(rows),
		band_(band),
		padding_(padding),
		scores_(TLS::get(scores_ptr)),
		hgap_(TLS::get(hgap_ptr))
	{
		scores_.clear();
		scores_.resize(2 * band + 1);
		hgap_.clear();
		hgap_.resize(2 * band + 2);
		hgap_front_ = &hgap_.front();
		score_front_ = &scores_.front();
	}

	inline void clear()
	{
		memset(hgap_front_, 0, (2 * band_ + 2) * sizeof(int));
		memset(score_front_, 0, (2 * band_ + 1) * sizeof(int));
	}

	inline Column_iterator begin(unsigned column)
	{
		if (column >= rows_ + padding_) {
			return Column_iterator(column, hgap_front_, score_front_, rows_ - band_, rows_, 0);
		}
		else if (column >= padding_) {
			unsigned pj(column - padding_);
			unsigned top_delta(pj >= band_ ? 0 : band_ - pj);
			unsigned query_start(pj >= band_ ? pj - band_ : 0);
			unsigned query_end(std::min(pj + band_ + 1, rows_));
			return Column_iterator(column, hgap_front_ + top_delta, score_front_ + top_delta, query_start, query_end, 1);
		}
		else {
			return Column_iterator(column, hgap_front_ + band_ + 1, score_front_ + band_ + 1, 0, band_, 0);
		}
	}

private:

	static TLS_PTR vector<int> *scores_ptr;
	static TLS_PTR vector<int> *hgap_ptr;

	const unsigned rows_, band_, padding_;
	int *hgap_front_, *score_front_;
	vector<int> &scores_, &hgap_;

};

TLS_PTR vector<int>* Padded_banded_DP_matrix::scores_ptr;
TLS_PTR vector<int>* Padded_banded_DP_matrix::hgap_ptr;

int smith_waterman(const sequence &query, const sequence &subject, unsigned band, unsigned padding, int op, int ep)
{
	unsigned qlen((unsigned)query.length());
	unsigned slen((unsigned)subject.length());
	Padded_banded_DP_matrix dp(slen, qlen, band, padding);

	int best = 0;
	dp.clear();

	for (unsigned j = 0; j < slen; ++j) {
		if (j < (unsigned)subject.clipping_offset_)
			continue;
		if (subject[j] == '\xff')
			break;
		Padded_banded_DP_matrix::Column_iterator it(dp.begin(j));
		int vgap = 0, hgap;
		while (!it.at_end()) {
			hgap = it.hgap();
			int current_cell = it.diag() + score_matrix(query[it.row_pos_], subject[j]);
			current_cell = std::max(current_cell, vgap);
			current_cell = std::max(current_cell, hgap);
			current_cell = std::max(current_cell, 0);
			best = std::max(best, current_cell);
			vgap -= ep;
			hgap -= ep;
			int open = current_cell - op;
			vgap = std::max(vgap, open);
			hgap = std::max(hgap, open);
			it.set_hgap(hgap);
			it.set_score(current_cell);
			++it;
		}
	}
	
	return best;
}
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

#include <vector>
#include "dp.h"
#include "../util/double_buffer.h"
#include "growing_buffer.h"
#include "../util/util.h"
#include "traceback.h"

using std::vector;
using std::pair;

template<typename _score>
struct Dp_matrix
{

	struct Column_iterator
	{

		inline Column_iterator(const pair<_score*, _score*> &score, _score* hgap, int query_len, int col) :
			score_(score),
			hgap_(hgap),
			end_(score_.second + query_len + 1),
			i_(0)
		{
			*score_.first = col == 0 ? 0 : -config.gap_open - col*config.gap_extend;
			++score_.second;
		}

		inline int row() const
		{
			return i_;
		}

		inline bool valid() const
		{
			return score_.second < end_;
		}

		inline _score& score()
		{
			return *score_.second;
		}

		inline _score diag() const
		{
			return *score_.first;
		}

		inline _score& hgap()
		{
			return *hgap_;
		}

		inline void operator++()
		{
			++i_;
			++score_.first;
			++score_.second;
			++hgap_;
		}

	private:
		pair<_score*, _score*> score_;
		_score* hgap_;
		const _score* const end_;
		int i_;

	};

	inline Column_iterator column(int j)
	{
		return Column_iterator(score_.get(), hgap_.data(), query_len_, j);
	}

	inline Dp_matrix(int query_len, int subject_len) :
		query_len_(query_len),
		score_(TLS::get(score_ptr)),
		hgap_(TLS::get(hgap_ptr))
	{
		score_.init(query_len + 1, subject_len + 1, 0);
		hgap_.clear();
		hgap_.insert(hgap_.end(), query_len, std::numeric_limits<int>::min() + 1);
		int *score = score_.last();
		int g = -config.gap_open - config.gap_extend;
		for (int i = 1; i <= query_len; ++i)
			score[i] = g--;
	}

	const Fixed_score_buffer<_score>& score_buffer() const
	{
		return score_;
	}

private:

	const int query_len_;
	Fixed_score_buffer<_score> &score_;
	vector<_score> &hgap_;
	static TLS_PTR Fixed_score_buffer<_score> *score_ptr;
	static TLS_PTR vector<_score> *hgap_ptr;

};

template<typename _score> TLS_PTR Fixed_score_buffer<_score>* Dp_matrix<_score>::score_ptr;
template<typename _score> TLS_PTR vector<_score>* Dp_matrix<_score>::hgap_ptr;

template<typename _score>
const Fixed_score_buffer<_score>& needleman_wunsch(sequence query, sequence subject, const _score& = int())
{
	using std::max;
	const int gap_open = config.gap_open + config.gap_extend, gap_extend = config.gap_extend;

	Dp_matrix<_score> mtx((unsigned)query.length(), (unsigned)subject.length());

	for (int j = 0; j < (int)subject.length(); ++j) {
		typename Dp_matrix<_score>::Column_iterator it = mtx.column(j);
		_score vgap = std::numeric_limits<int>::min() + 1;
		for (; it.valid(); ++it) {
			const _score match_score = score_matrix(subject[j], query[it.row()]);
			const _score s = max(max(it.diag() + match_score, vgap), it.hgap());
			const _score open = s - gap_open;
			vgap = max(vgap - gap_extend, open);
			it.hgap() = max(it.hgap() - gap_extend, open);
			it.score() = s;
		}
	}

	return mtx.score_buffer();
}

int needleman_wunsch(sequence query, sequence subject, int qbegin, int qend, int sbegin, int send, unsigned node, unsigned edge, vector<Diagonal_node> &diags, bool log)
{
	const sequence q = query.subseq(qbegin, qend), s = subject.subseq(sbegin, send);
	const Fixed_score_buffer<int> &dp = needleman_wunsch(q, s, int());
	Diagonal_node *d = &diags[node];
	unsigned start_node = d->edges[edge].node;
	Diagonal_node::Edge *f = &d->edges[edge];

	if (log)
		cout << dp << endl;

	const int gap_open = config.gap_open, gap_extend = config.gap_extend;
	int l, i = qend - qbegin, j = send - sbegin;
	const int score = dp(i, j);

	l = have_diag(dp, i, j, q, s, log);
	if (l > 0) {
		i -= l;
		j -= l;
		f->exact = true;
		f->j = sbegin + j;
	}

	while (i > 0 && j > 0) {
		if ((l = have_diag(dp, i, j, q, s, log)) > 0) {
			i -= l;
			j -= l;
			if (i != 0 || j != 0) {
				f->node = (unsigned)diags.size();
				diags.push_back(Diagonal_node(qbegin + i, sbegin + j, l, 0));
				f = &diags.back().edges[0];
				f->exact = true;
				f->j = sbegin + j;
			}
		}
		else if (have_hgap(dp, i, j, gap_open, gap_extend, l)) {
			j -= l;
		}
		else if (have_vgap(dp, i, j, gap_open, gap_extend, l)) {
			i -= l;
		}
		else
			throw std::runtime_error("Traceback error.");
	}

	f->node = start_node;
	return score;
}
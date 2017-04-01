/****
Copyright (c) 2016-2017, Benjamin Buchfink
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
#include "../output/output_format.h"

using std::vector;
using std::pair;

template<typename _score, typename _mode>
_score saturate(_score x)
{
	return x;
}

template<>
int saturate<int, Local>(int x)
{
	return std::max(x, 0);
}

template<typename _score,typename _mode>
void set_max_score(_score s, _score &max_score)
{
}

template<>
void set_max_score<int,Local>(int s, int &max_score)
{
	max_score = std::max(max_score, s);
}

template<typename _score, typename _mode>
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
			*score_.first = saturate<_score, _mode>(col == 0 ? 0 : -score_matrix.gap_open() - col*score_matrix.gap_extend());
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
		int g = -score_matrix.gap_open() - score_matrix.gap_extend();
		for (int i = 1; i <= query_len; ++i)
			score[i] = saturate<_score, _mode>(g--);
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

template<typename _score, typename _mode> TLS_PTR Fixed_score_buffer<_score>* Dp_matrix<_score,_mode>::score_ptr;
template<typename _score, typename _mode> TLS_PTR vector<_score>* Dp_matrix<_score,_mode>::hgap_ptr;

template<typename _score, typename _mode>
const Fixed_score_buffer<_score>& needleman_wunsch(sequence query, sequence subject, int &max_score, const _mode&, const _score&)
{
	using std::max;
	const int gap_open = score_matrix.gap_open() + score_matrix.gap_extend(), gap_extend = score_matrix.gap_extend();
	int m = 0;

	Dp_matrix<_score, _mode> mtx((unsigned)query.length(), (unsigned)subject.length());

	for (int j = 0; j < (int)subject.length(); ++j) {
		typename Dp_matrix<_score,_mode>::Column_iterator it = mtx.column(j);
		_score vgap = std::numeric_limits<int>::min() + 1;
		for (; it.valid(); ++it) {
			const _score match_score = score_matrix(subject[j], query[it.row()]);
			const _score s = saturate<_score, _mode>(max(max(it.diag() + match_score, vgap), it.hgap()));
			const _score open = s - gap_open;
			vgap = max(vgap - gap_extend, open);
			it.hgap() = max(it.hgap() - gap_extend, open);
			it.score() = s;
			set_max_score<_score, _mode>(s, m);
		}
	}

	max_score = m;
	return mtx.score_buffer();
}

int needleman_wunsch(sequence query, sequence subject, int qbegin, int qend, int sbegin, int send, unsigned node, unsigned edge, Diag_graph &diags, bool log)
{
	const sequence q = query.subseq(qbegin, qend), s = subject.subseq(sbegin, send);
	int max_score;
	const Fixed_score_buffer<int> &dp = needleman_wunsch(q, s, max_score, Global(), int());
	Diagonal_node *d = &diags[node];
	unsigned start_node = diags.edges[edge].node_out;
	vector<Diag_graph::Edge>::iterator f = diags.edges.begin() + edge;

	/*if (log)
		cout << dp << endl;*/

	const int gap_open = score_matrix.gap_open(), gap_extend = score_matrix.gap_extend();
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
				f->node_out = (unsigned)diags.nodes.size();
				diags.nodes.push_back(Diagonal_node(qbegin + i, sbegin + j, l, 0, diags.edges.size()));
				f = diags.add_edge(Diag_graph::Edge(0, sbegin + j, f->node_out, 0, true, Diagonal_node::finished, 0, 0));
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

	f->node_out = start_node;
	return score;
}

void smith_waterman(sequence q, sequence s, Hsp_data &out)
{
	int max_score;
	const Fixed_score_buffer<int> &dp = needleman_wunsch(q, s, max_score, Local(), int());
	pair<int, int> max_pos = dp.find(max_score);

	const int gap_open = score_matrix.gap_open(), gap_extend = score_matrix.gap_extend();
	int l, i = max_pos.first, j = max_pos.second, score;
	out.score = dp(i, j);
	out.query_range.end_ = i;
	out.subject_range.end_ = j;

	while ((score = dp(i, j)) > 0) {
		const int match_score = score_matrix(q[i - 1], s[j - 1]);
		if (score == match_score + dp(i - 1, j - 1)) {
			if (q[i - 1] == s[j - 1]) {
				out.transcript.push_back(op_match);
			}
			else {
				out.transcript.push_back(op_substitution, s[j - 1]);
			}
			--i;
			--j;
			++out.length;
		}
		else if (have_hgap(dp, i, j, gap_open, gap_extend, l)) {
			for (; l > 0; l--) {
				out.transcript.push_back(op_deletion, s[--j]);
				++out.length;
			}
		}
		else if (have_vgap(dp, i, j, gap_open, gap_extend, l)) {
			out.transcript.push_back(op_insertion, (unsigned)l);
			out.length += l;
			i -= l;
		}
		else
			throw std::runtime_error("Traceback error.");
	}

	out.query_range.begin_ = i;
	out.subject_range.begin_ = j;
	out.transcript.reverse();
	out.transcript.push_terminator();
}

void print_diag(int i0, int j0, int l, int score, const Diag_graph &diags, const sequence &query, const sequence &subject)
{
	Diagonal_segment ds(i0, j0, l, 0);
	unsigned n = 0;
	for (vector<Diagonal_node>::const_iterator d = diags.nodes.begin(); d != diags.nodes.end(); ++d) {
		if (d->intersect(ds).len > 0) {
			if (d->score == 0)
				continue;
			const int diff = score_range(query, subject, d->query_end(), d->subject_end(), j0 + l);
			if (n > 0)
				cout << "(";
			cout << "Diag n=" << d - diags.nodes.begin() << " i=" << i0 << " j=" << j0 << " len=" << l
				<< " prefix_score=" << score + score_range(query, subject, i0 + l, j0 + l, d->subject_end()) - std::min(diff, 0)
				<< " prefix_score2=" << diags.prefix_score(d - diags.nodes.begin(), j0 + l, 0);
			if (n > 0)
				cout << ")";
			cout << endl;
			++n;
		}
	}
	if(n == 0)
		cout << "Diag n=x i=" << i0 << " j=" << j0 << " len=" << l << " prefix_score=" << score << endl;
}

void smith_waterman(sequence q, sequence s, const Diag_graph &diags, Text_buffer &buf)
{
	Hsp_data hsp;
	smith_waterman(q, s, hsp);
	Hsp_data::Iterator i = hsp.begin();
	int i0 = -1, j0 = -1, l = 0, score = 0;
	for (; i.good(); ++i) {
		switch (i.op()) {
		case op_match:
		case op_substitution:
			if (i0 < 0) {
				i0 = i.query_pos;
				j0 = i.subject_pos;
				l = 0;
			}
			score += score_matrix(q[i.query_pos], s[i.subject_pos]);
			++l;
			break;
		case op_deletion:
		case op_insertion:
			if (i0 >= 0) {
				print_diag(i0, j0, l, score, diags, q, s);
				score -= score_matrix.gap_open() + score_matrix.gap_extend();
				i0 = -1;
				j0 = -1;
			}
			else
				score -= score_matrix.gap_extend();
		}
	}
	print_diag(i0, j0, l, score, diags, q, s);
	buf.clear();
	Pairwise_format().print_match(Hsp_context(hsp, 0, q, q, "", 0, 0, "", 0, 0, 0), buf);
}
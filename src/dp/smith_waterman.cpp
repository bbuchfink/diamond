/****
Copyright (c) 2016, University of Tuebingen, Benjamin Buchfink
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

using std::vector;
using std::pair;

template<typename _score>
struct Dp_matrix
{

	struct Column_iterator
	{

		inline Column_iterator(const pair<_score*, _score*> &score, const pair<_score*, _score*> &hgap, int query_len):
			score_(score),
			hgap_(hgap),
			end_(score_.second + query_len),
			i_(0)
		{
			*score_.first = 0;
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

		inline _score hgap_in() const
		{
			return *hgap_.first;
		}

		inline _score& hgap_out()
		{
			return *hgap_.second;
		}

		inline void operator++()
		{
			++i_;
			++score_.first;
			++score_.second;
			++hgap_.first;
			++hgap_.second;
		}

	private:
		pair<_score*, _score*> score_, hgap_;
		const _score* const end_;
		int i_;

	};

	inline Column_iterator column(int j)
	{
		return Column_iterator(score_.get(), hgap_.get(int()), query_len_);
	}

	inline Dp_matrix(int query_len, int subject_len):
		query_len_(query_len),
		score_(TLS::get(score_ptr)),
		hgap_(TLS::get(hgap_ptr))
	{
		score_.init(query_len+1, subject_len+1, 0);
		hgap_.init(query_len, 0, 0, 0);
	}

	const Fixed_score_buffer<_score>& score_buffer() const
	{
		return score_;
	}
	
private:

	const int query_len_;
	Fixed_score_buffer<_score> &score_;
	Double_buffer<_score> &hgap_;
	static TLS_PTR Fixed_score_buffer<_score> *score_ptr;
	static TLS_PTR Double_buffer<_score> *hgap_ptr;

};

template<typename _score>
void smith_waterman(const Letter *query, unsigned query_len, local_match &hssp, _score gap_open, _score gap_extend, vector<char> &transcript_buf, const _score& = int())
{
	using std::max;

	_score max_score = 0, column_max = 0;
	int j = 0, i_max = -1, j_best = -1, i_best = -1;
	Dp_matrix<_score> mtx(query_len, hssp.total_subject_len_);
	const Letter *subject = hssp.subject_;

	while (*subject != '\xff') {
		typename Dp_matrix<_score>::Column_iterator it = mtx.column(j);
		_score vgap = 0;		
		for (; it.valid(); ++it) {
			const _score match_score = score_matrix(mask_critical(*subject), query[it.row()]);
			const _score s = max(max(it.diag() + match_score, vgap), it.hgap_in());
			s = max(s, 0);
			if (s > column_max) {
				column_max = s;
				i_max = it.row();
			}
			const _score open = s - gap_open;
			vgap = max(vgap - gap_extend, open);
			it.hgap_out() = max(it.hgap_in() - gap_extend, open);
			it.score() = s;
		}

		if (column_max > max_score) {
			max_score = column_max;
			j_best = j;
			i_best = i_max;
		}
		++subject;
		++j;
	}

	return traceback<_dir, _score>(query, subject, mtx.score_buffer(), band, gap_open, gap_extend, j_best, i_best, max_score, transcript_buf);
}

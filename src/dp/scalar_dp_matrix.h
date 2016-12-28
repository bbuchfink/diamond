/****
Copyright (c) 2014, University of Tuebingen
Author: Benjamin Buchfink
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

#ifndef SCALAR_DP_MATRIX_H_
#define SCALAR_DP_MATRIX_H_

#include <vector>
#include "../util/double_buffer.h"
#include "growing_buffer.h"
#include "../util/util.h"
#include "floating_sw.h"

using std::vector;
using std::pair;

template<typename _score, typename _traceback>
struct Score_buffer { };

template<typename _score>
struct Score_buffer<_score,Score_only>
{
	typedef Double_buffer<_score> Type;
};

template<typename _score>
struct Score_buffer<_score,Traceback>
{
	typedef Growing_buffer<_score> Type;
};

template<typename _score, typename _traceback>
struct Scalar_dp_matrix
{

	struct Column_iterator
	{

		inline Column_iterator(const pair<_score*,_score*> &score, const pair<_score*,_score*> &hgap, int j, int i, int delta, int band):
			score_ (score),
			hgap_ (hgap),
			end_ (score_.second + 2*band + 1),
			i_ (std::max(i - band, 0))
		{
			assert(delta >= 0 && j >= 0 && i >= 0 && band >= 0);
			if(j == 0)
				score_.first[band] = 0;
			const int offset = i_ - i + band;
			score_.second += offset;
			hgap_.second += offset;
			hgap_.first += offset + delta;
			score_.first += offset + delta - 1;
		}

		inline int row() const
		{ return i_; }

		inline bool valid() const
		{ return score_.second < end_; }

		inline _score& score()
		{ return *score_.second; }

		inline _score diag() const
		{ return *score_.first; }

		inline _score hgap_in() const
		{ return *hgap_.first; }

		inline _score& hgap_out()
		{ return *hgap_.second; }

		inline void operator++()
		{
			++i_;
			++score_.first;
			++score_.second;
			++hgap_.first;
			++hgap_.second;
		}

	private:
		pair<_score*,_score*> score_, hgap_;
		const _score* const end_;
		int i_;

	};

	inline Column_iterator column(int j, int i_max)
	{
		int i = std::max(current_i_, i_max+1), delta = i - current_i_;
		current_i_ = i;
		return Column_iterator (score_.get(i), hgap_.get(int ()), j, i, delta, band_);
	}

	inline Scalar_dp_matrix(int band):
		band_ (band),
		band_max_ (2*band+1),
		current_i_ (-1),
		score_ (TLS::get(score_ptr)),
		hgap_ (TLS::get(hgap_ptr))
	{
		score_.init(band_max_, band_+1, 1, minus_inf);
		hgap_.init(band_max_, band_+1, 1, minus_inf);
	}

	const typename Score_buffer<_score,_traceback>::Type& score_buffer() const
	{ return score_; }

	enum {
		minus_inf = -65536
	};

private:

	const int band_, band_max_;
	int current_i_;
	typename Score_buffer<_score,_traceback>::Type &score_;
	Double_buffer<_score> &hgap_;
	static TLS_PTR typename Score_buffer<_score,_traceback>::Type *score_ptr;
	static TLS_PTR Double_buffer<_score> *hgap_ptr;

};

#endif /* SCALAR_DP_MATRIX_H_ */

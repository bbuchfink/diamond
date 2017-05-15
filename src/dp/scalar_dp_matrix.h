/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
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

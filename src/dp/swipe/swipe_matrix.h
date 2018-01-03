/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

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

#ifndef SWIPE_MATRIX_H_
#define SWIPE_MATRIX_H_

#include <vector>
#include <string.h>
#include "../score_vector.h"
#include "../../util/tls.h"

template<typename _score>
struct SwipeMatrix
{
	typedef score_vector<_score> sv;
	struct ColumnIterator
	{
		ColumnIterator(sv* hgap_front, sv* score_front) :
			hgap_ptr_(hgap_front),
			score_ptr_(score_front)
		{ }
		inline void operator++()
		{
			++hgap_ptr_; ++score_ptr_;
		}
		inline sv hgap() const
		{
			return *hgap_ptr_;
		}
		inline sv diag() const
		{
			return *score_ptr_;
		}
		inline void set_hgap(const sv& x)
		{
			*hgap_ptr_ = x;
		}
		inline void set_score(const sv& x)
		{
			*score_ptr_ = x;
		}
		sv *hgap_ptr_, *score_ptr_;
	};
	SwipeMatrix(int rows) :
		hgap_(TLS::get(hgap_ptr)),
		score_(TLS::get(score_ptr))
	{
		hgap_.clear();
		hgap_.resize(rows);
		score_.clear();
		score_.resize(rows + 1);
		memset(hgap_.data(), 0, rows * sizeof(sv));
		memset(score_.data(), 0, (rows + 1) * sizeof(sv));
	}
	inline ColumnIterator begin()
	{
		return ColumnIterator(&hgap_[0], &score_[0]);
	}
	void set_zero(int c)
	{
		const int l = (int)hgap_.size();
		for (int i = 0; i < l; ++i) {
			hgap_[i].set(c, 0);
			score_[i].set(c, 0);
		}
		score_[l].set(c, 0);
	}
private:
	std::vector<sv> &hgap_, &score_;
	static TLS_PTR std::vector<sv> *hgap_ptr, *score_ptr;
};


template<typename _score>
struct BandedSwipeMatrix
{
	typedef score_vector<_score> sv;
	struct ColumnIterator
	{
		ColumnIterator(sv* hgap_front, sv* score_front) :
			hgap_ptr_(hgap_front),
			score_ptr_(score_front)
		{ }
		inline void operator++()
		{
			++hgap_ptr_; ++score_ptr_;
		}
		inline sv hgap() const
		{
			return *(hgap_ptr_ + 1);
		}
		inline sv diag() const
		{
			return *score_ptr_;
		}
		inline void set_hgap(const sv& x)
		{
			*hgap_ptr_ = x;
		}
		inline void set_score(const sv& x)
		{
			*score_ptr_ = x;
		}
		sv *hgap_ptr_, *score_ptr_;
	};
	BandedSwipeMatrix(int band) :
		hgap_(TLS::get(hgap_ptr)),
		score_(TLS::get(score_ptr))
	{
		hgap_.clear();
		hgap_.resize(band + 1);
		score_.clear();
		score_.resize(band);
		memset(hgap_.data(), 0, (band + 1) * sizeof(sv));
		memset(score_.data(), 0, rows * sizeof(sv));
	}
	inline ColumnIterator begin(int offset)
	{
		return ColumnIterator(&hgap_[offset], &score_[offset]);
	}
private:
	std::vector<sv> &hgap_, &score_;
	static TLS_PTR std::vector<sv> *hgap_ptr, *score_ptr;
};

#endif
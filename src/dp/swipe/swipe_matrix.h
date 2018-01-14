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
#include <stdlib.h>
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
		/*memset(hgap_.data(), 0, (band + 1) * sizeof(sv));
		memset(score_.data(), 0, band * sizeof(sv));*/
		size_t i = 0;
		sv z = sv();
		for (; i < band; ++i) {
			hgap_[i] = z;
			score_[i] = z;
		}
		hgap_[i] = z;
	}
	inline ColumnIterator begin(int offset)
	{
		return ColumnIterator(&hgap_[offset], &score_[offset]);
	}
private:
	std::vector<sv> &hgap_, &score_;
	static TLS_PTR std::vector<sv> *hgap_ptr, *score_ptr;
};

template<typename _score>
struct BandedSwipeTracebackMatrix
{
	typedef score_vector<_score> sv;
	struct ColumnIterator
	{
		ColumnIterator(sv* hgap_front, sv* score_front, sv* hgap_front1, sv* score_front1) :
			hgap_ptr_(hgap_front),
			score_ptr_(score_front),
			hgap_ptr1_(hgap_front1),
			score_ptr1_(score_front1)
		{ }
		inline void operator++()
		{
			++hgap_ptr_; ++score_ptr_; ++hgap_ptr1_; ++score_ptr1_;
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
			*hgap_ptr1_ = x;
		}
		inline void set_score(const sv& x)
		{
			*score_ptr1_ = x;
		}
		void set_zero()
		{
			(score_ptr1_ - 1)->zero();
		}
		sv *hgap_ptr_, *score_ptr_, *hgap_ptr1_, *score_ptr1_;
	};
	BandedSwipeTracebackMatrix(size_t band, size_t cols) :
		band_ (band),
		hgap_(TLS::get(hgap_ptr)),
		score_(TLS::get(score_ptr))
	{
		hgap_.clear();
		hgap_.resize((band + 1) * (cols + 1));
		score_.clear();
		score_.resize(band * (cols + 1));
		size_t i = 0;
		sv z;
		z.zero();
		for (; i < band; ++i) {
			hgap_[i] = z;
			score_[i] = z;
		}
		for (i = 0; i < cols; ++i)
			hgap_[i*(band + 1) + band] = z;
	}
	inline ColumnIterator begin(size_t offset, size_t col)
	{
		return ColumnIterator(&hgap_[col*(band_ + 1) + offset + 1], &score_[col*band_ + offset], &hgap_[(col + 1)*(band_ + 1) + offset], &score_[(col + 1)*band_ + offset]);
	}
private:
	const size_t band_;
	std::vector<sv> &hgap_, &score_;
	static TLS_PTR std::vector<sv> *hgap_ptr, *score_ptr;
};

template<typename _score>
struct Banded3FrameSwipeTracebackMatrix
{

	typedef score_vector<_score> sv;

	struct ColumnIterator
	{
		ColumnIterator(sv* hgap_front, sv* score_front, sv* score_front1) :
			hgap_ptr_(hgap_front),
			score_ptr_(score_front),
			score_ptr1_(score_front1)
		{
			sm4.zero();
			sm3 = *(score_ptr_++);
			sm2 = *(score_ptr_);
		}
		inline void operator++()
		{
			++hgap_ptr_; ++score_ptr_; ++score_ptr1_;
			sm4 = sm3;
			sm3 = sm2;
			sm2 = *score_ptr_;
		}
		inline sv hgap() const
		{
			return *(hgap_ptr_ + 3);
		}
		inline void set_hgap(const sv& x)
		{
			*hgap_ptr_ = x;
		}
		inline void set_score(const sv& x)
		{
			*score_ptr1_ = x;
		}
		void set_zero()
		{
			(score_ptr1_ - 1)->zero();
			(score_ptr1_ - 2)->zero();
			(score_ptr1_ - 3)->zero();
		}
		sv *hgap_ptr_, *score_ptr_, *score_ptr1_;
		sv sm4, sm3, sm2;
	};

	struct TracebackIterator
	{
		TracebackIterator(const _score *score, size_t band, int frame, int i, int j):
			band_(band),
			score_(score),
			frame(frame),
			i(i),
			j(j)
		{}
		_score score() const
		{
			return *score_;
		}
		_score sm3() const
		{
			return *(score_ - (band_ + 1) * 8);
		}
		_score sm4() const
		{
			return *(score_ - (band_ + 2) * 8);
		}
		_score sm2() const
		{
			return *(score_ - band_ * 8);
		}
		void walk_diagonal()
		{
			score_ -= (band_ + 1) * 8;
			--i;
			--j;
			assert(i >= -1 && j >= -1);
		}
		void walk_forward_shift()
		{
			score_ -= (band_ + 2) * 8;
			--i;
			--j;
			--frame;
			if (frame == -1) {
				frame = 2;
				--i;
			}
			assert(i >= -1 && j >= -1);
		}
		void walk_reverse_shift()
		{
			score_ -= band_ * 8;
			--i;
			--j;
			++frame;
			if (frame == 3) {
				frame = 0;
				++i;
			}
			assert(i >= -1 && j >= -1);
		}
		int walk_gap(int d0, int d1)
		{
			const int i0 = std::max(d0 + j, 0), j0 = std::max(i - d1, -1);
			const _score *h = score_ - (band_ - 2) * 8, *h0 = score_ - (j - j0) * (band_ - 2) * 8;
			const _score *v = score_ - 3 * 8, *v0 = score_ - (i - i0 + 1) * 3 * 8;
			const _score score = this->score();
			const _score e = score_matrix.gap_extend();
			_score g = score_matrix.gap_open() + e;
			int l = 1;
			while (v>v0 && h>h0) {
				if (score + g == *h) {
					walk_hgap(h, l);
					return -l;
				}
				else if (score + g == *v) {
					walk_vgap(v, l);
					return l;
				}
				h -= (band_ - 2) * 8;
				v -= 3 * 8;
				++l;
				g += e;
			}
			while (v>v0) {
				if (score + g == *v) {
					walk_vgap(v, l);
					return l;
				}
				v -= 3 * 8;
				++l;
				g += e;
			}
			while (h>h0) {
				if (score + g == *h) {
					walk_hgap(h, l);
					return -l;
				}
				h -= (band_ - 2) * 8;
				++l;
				g += e;
			}
			throw std::runtime_error("Traceback error.");
		}
		void walk_hgap(const _score *h, int l)
		{
			score_ = h;
			j -= l;
			assert(i >= -1 && j >= -1);
		}
		void walk_vgap(const _score *v, int l)
		{
			score_ = v;
			i -= l;
			assert(i >= -1 && j >= -1);
		}
		const size_t band_;
		const _score *score_;
		int frame, i, j;
	};

	TracebackIterator traceback(size_t col, int i0, int j, size_t channel, _score score) const
	{
		const _score *s = (_score*)(&score_[col*(band_ + 1)]) + channel;
		// TODO: Clipping
		for (int i = 0; i < band_; ++i, s += 8)
			if (*s == score)
				return TracebackIterator(s, band_, i % 3, i0 + i / 3, j);
		throw std::runtime_error("Trackback error.");
	}

	Banded3FrameSwipeTracebackMatrix(size_t band, size_t cols) :
		band_(band),
		hgap_(TLS::get(hgap_ptr)),
		score_(TLS::get(score_ptr))
	{
		hgap_.clear();
		hgap_.resize(band + 3);
		score_.clear();
		score_.resize((band + 1) * (cols + 1));
		size_t i = 0;
		sv z;
		z.zero();
		for (; i < band + 1; ++i) {
			hgap_[i] = z;
			score_[i] = z;
		}
		hgap_[i++] = z;
		hgap_[i] = z;
		for (i = 0; i < cols; ++i) {
			score_[i*(band + 1) + band] = z;
		}
	}

	inline ColumnIterator begin(size_t offset, size_t col)
	{
		return ColumnIterator(&hgap_[offset], &score_[col*(band_ + 1) + offset], &score_[(col + 1)*(band_ + 1) + offset]);
	}

private:
	const size_t band_;
	std::vector<sv> &hgap_, &score_;
	static TLS_PTR std::vector<sv> *hgap_ptr, *score_ptr;

};

#endif
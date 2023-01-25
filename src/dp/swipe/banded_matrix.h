/****
DIAMOND protein aligner
Copyright (C) 2016-2021 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
						
Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#include <utility>
#include "../../util/data_structures/mem_buffer.h"

using std::pair;

namespace DP { namespace BandedSwipe {
namespace DISPATCH_ARCH {

template<typename Sv>
struct Matrix
{
	using Score = typename ::DISPATCH_ARCH::ScoreTraits<Sv>::Score;
	static constexpr int CHANNELS = ::DISPATCH_ARCH::ScoreTraits<Sv>::CHANNELS;
	struct ColumnIterator
	{
		ColumnIterator(Sv* hgap_front, Sv* score_front) :
			hgap_ptr_(hgap_front),
			score_ptr_(score_front)
		{ }
		inline void operator++()
		{
			++hgap_ptr_; ++score_ptr_;
		}
		inline Sv hgap() const
		{
			return *(hgap_ptr_ + 1);
		}
		inline Sv diag() const
		{
			return *score_ptr_;
		}
		inline void set_hgap(const Sv& x)
		{
			*hgap_ptr_ = x;
		}
		inline void set_score(const Sv& x)
		{
			*score_ptr_ = x;
		}
		std::nullptr_t stat() {
			return nullptr;
		}
		std::nullptr_t hstat() {
			return nullptr;
		}
		std::nullptr_t trace_mask() {
			return nullptr;
		}
		void set_hstat(std::nullptr_t) {}
		inline void set_zero() {}
		Sv *hgap_ptr_, *score_ptr_;
	};
	Matrix(int band, size_t cols, Sv init = Sv()) :
		band_(band)
	{
		hgap_.resize(band + 1);
		score_.resize(band);
		std::fill(hgap_.begin(), hgap_.end(), init);
		std::fill(score_.begin(), score_.end(), init);

	}
	void init_channel_diag(int channel, int offset) {
		Score* ptr = (Score*)score_.begin();
		ptr[offset * CHANNELS + channel] = 0;
	}
	void init_channel_nw(int channel, int offset, ::Score gap_open, ::Score gap_extend) {
		Score* ptr = (Score*)score_.begin();
		ptr[offset * CHANNELS + channel] = 0;
		//Score s = -gap_open;
		for (int i = offset - 1; i >= 0; --i) {
			//s -= gap_extend;
			//ptr[i * CHANNELS + channel] = s;
			ptr[i * CHANNELS + channel] = std::numeric_limits<Score>::min();
		}
		//s = -gap_open;
		for (int i = offset + 1; i < (int)score_.size(); ++i) {
			//s -= gap_extend;
			//ptr[i * CHANNELS + channel] = s;
			ptr[i * CHANNELS + channel] = std::numeric_limits<Score>::min();
		}
		ptr = (Score*)hgap_.begin();
		for (size_t i = 0; i < hgap_.size(); ++i)
			ptr[i * CHANNELS + channel] = std::numeric_limits<Score>::min();
	}
	void init_channels_nw(int offset, ::Score gap_open, ::Score gap_extend) {
		Score* ptr = (Score*)score_.begin();
		std::fill(&ptr[offset * CHANNELS], &ptr[(offset+1) * CHANNELS], 0);
		Score s = -gap_open;
		for (int i = offset - 1; i >= 0; --i) {
			s -= gap_extend;
			std::fill(&ptr[i * CHANNELS], &ptr[(i + 1) * CHANNELS], s);
		}
		s = -gap_open;
		for (int i = offset + 1; i < (int)score_.size(); ++i) {
			s -= gap_extend;
			std::fill(&ptr[i * CHANNELS], &ptr[(i + 1) * CHANNELS], s);
		}
		ptr = (Score*)hgap_.begin();
		std::fill(ptr, ptr + hgap_.size() * CHANNELS, std::numeric_limits<Score>::min());
	}
	inline ColumnIterator begin(int offset, int col)
	{
		return ColumnIterator(&hgap_[offset], &score_[offset]);
	}
	int band() const {
		return band_;
	}
	Sv operator[](int i) const {
		return score_[i];
	}
#ifdef __APPLE__
	MemBuffer<Sv> hgap_, score_;
#else
	static thread_local MemBuffer<Sv> hgap_, score_;
#endif
private:
	int band_;	
};

template<typename _sv>
struct TracebackMatrix
{

	typedef void* Stat;
	typedef typename ::DISPATCH_ARCH::ScoreTraits<_sv>::Score Score;
	static constexpr int CHANNELS = ::DISPATCH_ARCH::ScoreTraits<_sv>::CHANNELS;

	struct ColumnIterator
	{
		ColumnIterator(_sv* hgap_front, _sv* score_front, _sv* score_front1) :
			hgap_ptr_(hgap_front),
			score_ptr_(score_front),
			score_ptr1_(score_front1)
		{ }
		inline void operator++()
		{
			++hgap_ptr_; ++score_ptr_; ++score_ptr1_;
		}
		inline _sv hgap() const
		{
			return *(hgap_ptr_ + 1);
		}
		inline _sv diag() const
		{
			return *score_ptr_;
		}
		inline void set_hgap(const _sv& x)
		{
			*hgap_ptr_ = x;
		}
		inline void set_score(const _sv& x)
		{
			*score_ptr1_ = x;
		}
		void set_zero()
		{
			*(score_ptr1_ - 1) = _sv();
		}
		std::nullptr_t stat() {
			return nullptr;
		}
		std::nullptr_t hstat() {
			return nullptr;
		}
		std::nullptr_t trace_mask() {
			return nullptr;
		}
		void set_hstat(std::nullptr_t) {}
		_sv *hgap_ptr_, *score_ptr_, *score_ptr1_;
	};

	struct TracebackIterator
	{
		TracebackIterator(const Score *score, size_t band, int i, int j) :
			band_(band),
			score_(score),
			i(i),
			j(j)
		{
			assert(i >= 0 && j >= 0);
		}
		Score score() const
		{
			return *score_;
		}
		Score diag() const
		{
			return *(score_ - band_ * CHANNELS);
		}
		void walk_diagonal()
		{
			score_ -= band_ * ::DISPATCH_ARCH::ScoreTraits<_sv>::CHANNELS;
			--i;
			--j;
			assert(i >= -1 && j >= -1);
		}
		pair<Edit_operation, int> walk_gap(int d0, int d1)
		{
			const int i0 = std::max(d0 + j, 0), j0 = std::max(i - d1, -1);
			const Score *h = score_ - (band_ - 1) * CHANNELS, *h0 = score_ - (j - j0) * (band_ - 1) * CHANNELS;
			const Score *v = score_ - CHANNELS, *v0 = score_ - (i - i0 + 1) * CHANNELS;
			const Score e = score_matrix.gap_extend();
			Score score = this->score() + (Score)score_matrix.gap_open() + e;
			int l = 1;
			while (v > v0 && h > h0) {
				if (score == *h) {
					walk_hgap(h, l);
					return std::make_pair(op_deletion, l);
				}
				else if (score == *v) {
					walk_vgap(v, l);
					return std::make_pair(op_insertion, l);
				}
				h -= (band_ - 1) * CHANNELS;
				v -= CHANNELS;
				++l;
				score += e;
			}
			while (v > v0) {
				if (score == *v) {
					walk_vgap(v, l);
					return std::make_pair(op_insertion, l);
				}
				v -= CHANNELS;
				++l;
				score += e;
			}
			while (h > h0) {
				if (score == *h) {
					walk_hgap(h, l);
					return std::make_pair(op_deletion, l);
				}
				h -= (band_ - 1) * CHANNELS;
				++l;
				score += e;
			}
			throw std::runtime_error("Traceback error.");
		}
		void walk_hgap(const Score *h, int l)
		{
			score_ = h;
			j -= l;
			assert(i >= -1 && j >= -1);
		}
		void walk_vgap(const Score *v, int l)
		{
			score_ = v;
			i -= l;
			assert(i >= -1 && j >= -1);
		}
		const size_t band_;
		const Score *score_;
		int i, j;
	};

	TracebackIterator traceback(size_t col, int i0, int j, int query_len, size_t channel, Score score) const
	{
		const int i_ = std::max(-i0, 0),
			i1 = (int)std::min(band_, size_t(query_len - i0));
		const Score *s = (Score*)(&score_[col*band_ + i_]) + channel;
		for (int i = i_; i < i1; ++i, s += CHANNELS)
			if (*s == score)
				return TracebackIterator(s, band_, i0 + i, j);
		throw std::runtime_error("Trackback error.");
	}

	TracebackMatrix(size_t band, size_t cols) :
		band_(band)
	{
		hgap_.resize(band + 1);
		score_.resize(band * (cols + 1));
		std::fill(hgap_.begin(), hgap_.end(), _sv());
		std::fill(score_.begin(), score_.begin() + band, _sv());
	}

	inline ColumnIterator begin(size_t offset, size_t col)
	{
		return ColumnIterator(&hgap_[offset], &score_[col*band_ + offset], &score_[(col + 1)*band_ + offset]);
	}

	_sv operator[](int i) const {
		return _sv();
	}

	MemBuffer<_sv> hgap_, score_;

private:

	const size_t band_;

};

template<typename _sv>
struct TracebackVectorMatrix
{
	typedef typename ::DISPATCH_ARCH::ScoreTraits<_sv>::TraceMask TraceMask;
	typedef void* Stat;
	struct ColumnIterator
	{
		ColumnIterator(_sv* hgap_front, _sv* score_front, TraceMask* trace_mask_front) :
			hgap_ptr_(hgap_front),
			score_ptr_(score_front),
			trace_mask_ptr_(trace_mask_front)
		{ }
		inline void operator++()
		{
			++hgap_ptr_; ++score_ptr_; ++trace_mask_ptr_;
		}
		inline _sv hgap() const
		{
			return *(hgap_ptr_ + 1);
		}
		inline _sv diag() const
		{
			return *score_ptr_;
		}
		inline TraceMask* trace_mask() {
			return trace_mask_ptr_;
		}
		inline void set_hgap(const _sv& x)
		{
			*hgap_ptr_ = x;
		}
		inline void set_score(const _sv& x)
		{
			*score_ptr_ = x;
		}
		std::nullptr_t stat() {
			return nullptr;
		}
		std::nullptr_t hstat() {
			return nullptr;
		}
		void set_hstat(std::nullptr_t) {}
		inline void set_zero() {}
		_sv *hgap_ptr_, *score_ptr_;
		TraceMask* trace_mask_ptr_;
	};

	struct TracebackIterator
	{
		TracebackIterator(const TraceMask *mask, size_t band, int i, int j, const int channel) :
			band_(band),
			mask_(mask),
			channel_mask_vgap(TraceMask::vmask(channel)),
			channel_mask_hgap(TraceMask::hmask(channel)),
			i(i),
			j(j)
		{
			assert(i >= 0 && j >= 0);
		}
		TraceMask mask() const {
			return *mask_;
		}
		void walk_diagonal()
		{
			mask_ -= band_;
			--i;
			--j;
			assert(i >= -1 && j >= -1);
		}
		pair<Edit_operation, int> walk_gap()
		{
			if (mask_->gap & channel_mask_vgap) {
				int l = 0;
				do {
					++l;
					--i;
					--mask_;
				} while (((mask_->open & channel_mask_vgap) == 0) && (i > 0));
				return std::make_pair(op_insertion, l);
			}
			else {
				int l = 0;
				do {
					++l;
					--j;
					mask_ -= band_ - 1;
				} while (((mask_->open & channel_mask_hgap) == 0) && (j > 0));
				return std::make_pair(op_deletion, l);
			}
		}
		const size_t band_;
		const TraceMask* mask_;
		const decltype(TraceMask::gap) channel_mask_vgap, channel_mask_hgap;
		int i, j;
	};

	TracebackIterator traceback(size_t col, int i0, int band_i, int j, int query_len, const int channel) const
	{
		return TracebackIterator(&trace_mask_[col*band_ + band_i], band_, i0 + band_i, j, channel);
	}

	TracebackVectorMatrix(int band, size_t cols) :
		band_(band)
	{
		hgap_.resize(band + 1);
		score_.resize(band);
		trace_mask_.resize((cols + 1) * band);
		std::fill(hgap_.begin(), hgap_.end(), _sv());
		std::fill(score_.begin(), score_.end(), _sv());
	}
	
	inline ColumnIterator begin(int offset, int col)
	{
		return ColumnIterator(&hgap_[offset], &score_[offset], &trace_mask_[size_t(col + 1)*(size_t)band_ + (size_t)offset]);
	}

	int band() const {
		return band_;
	}

	_sv operator[](int i) const {
		return _sv();
	}

#ifdef __APPLE__
	MemBuffer<_sv> hgap_, score_;
#else
	static thread_local MemBuffer<_sv> hgap_, score_;
#endif
	MemBuffer<TraceMask> trace_mask_;
private:
	int band_;
};

#ifndef __APPLE__
template<typename Sv> thread_local MemBuffer<Sv> Matrix<Sv>::hgap_;
template<typename Sv> thread_local MemBuffer<Sv> Matrix<Sv>::score_;
template<typename Sv> thread_local MemBuffer<Sv> TracebackVectorMatrix<Sv>::hgap_;
template<typename Sv> thread_local MemBuffer<Sv> TracebackVectorMatrix<Sv>::score_;
#endif

template<typename Sv, bool Traceback>
struct SelectMatrix {
};

template<typename Sv>
struct SelectMatrix<Sv, true> {
	using Type = TracebackVectorMatrix<Sv>;
};

template<typename Sv>
struct SelectMatrix<Sv, false> {
	using Type = Matrix<Sv>;
};

}}}
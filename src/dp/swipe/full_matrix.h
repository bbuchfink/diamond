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

namespace DP { namespace Swipe {
namespace DISPATCH_ARCH {

template<typename _sv>
struct Matrix
{
	struct ColumnIterator
	{
		ColumnIterator(_sv* hgap_front, _sv* score_front) :
			hgap_ptr_(hgap_front),
			score_ptr_(score_front)
		{ }
		inline void operator++()
		{
			++hgap_ptr_; ++score_ptr_;
		}
		inline _sv hgap() const
		{
			return *hgap_ptr_;
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
			*score_ptr_ = x;
		}
		std::nullptr_t trace_mask() {
			return nullptr;
		}
		_sv *hgap_ptr_, *score_ptr_;
	};
	Matrix(int rows, int)
	{
		hgap_.resize(rows);
		score_.resize(rows + 1);
		const auto z = _sv();
		std::fill(hgap_.begin(), hgap_.end(), z);
		std::fill(score_.begin(), score_.end(), z);
	}
	inline ColumnIterator begin(int)
	{
		return ColumnIterator(hgap_.begin(), score_.begin());
	}
	void set_zero(int c)
	{
		const int l = (int)hgap_.size();
		const auto z = extract_channel(_sv(), 0);
		for (int i = 0; i < l; ++i) {
			set_channel(hgap_[i], c, z);
			set_channel(score_[i], c, z);
		}
		set_channel(score_[l], c, z);
	}
	constexpr int cols() const {
		return 1;
	}
	_sv operator[](int i) const {
		return score_[i + 1];
	}
private:
#ifdef __APPLE__
	MemBuffer<_sv> hgap_, score_;
#else
	static thread_local MemBuffer<_sv> hgap_, score_;
#endif
};

#ifndef __APPLE__
template<typename _sv> thread_local MemBuffer<_sv> Matrix<_sv>::hgap_;
template<typename _sv> thread_local MemBuffer<_sv> Matrix<_sv>::score_;
#endif

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
			return *hgap_ptr_;
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
		TracebackIterator(const TraceMask *mask, const TraceMask* mask_begin, const TraceMask* mask_end, int rows, int i, int j, int channel) :
			rows_(rows),
			mask_(mask),
			mask_begin_(mask_begin),
			mask_end_(mask_end),
			channel_mask_vgap(TraceMask::vmask(channel)),
			channel_mask_hgap(TraceMask::hmask(channel)),
			i(i),
			j(j)
		{
			assert(i >= 0 && j >= 0);
		}
		void wrap_mask() {
			if (mask_ < mask_begin_)
				mask_ = mask_end_ - (mask_begin_ - mask_);
		}
		TraceMask mask() const {
			return *mask_;
		}
		void walk_diagonal()
		{
			mask_ -= rows_ + 1;
			wrap_mask();
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
					mask_ -= rows_;
					wrap_mask();
				} while (((mask_->open & channel_mask_hgap) == 0) && (j > 0));
				return std::make_pair(op_deletion, l);
			}
		}
		const int rows_;
		const TraceMask* mask_, *mask_begin_, *mask_end_;
		const decltype(TraceMask::gap) channel_mask_vgap, channel_mask_hgap;
		int i, j;
	};

	TracebackIterator traceback(int col, int i, int j, int channel) const
	{
		return TracebackIterator(&trace_mask_[col*rows_ + i], trace_mask_.begin(), trace_mask_.end(), rows_, i, j, channel);
	}

	TracebackVectorMatrix(int rows, int cols) :
		rows_(rows),
		cols_(cols)
	{
		hgap_.resize(rows);
		score_.resize(rows + 1);
		trace_mask_.resize(cols * rows);
		std::fill(hgap_.begin(), hgap_.end(), _sv());
		std::fill(score_.begin(), score_.end(), _sv());
	}

	inline ColumnIterator begin(int col)
	{
		return ColumnIterator(hgap_.begin(), score_.begin(), &trace_mask_[col*rows_]);
	}

	void set_zero(int c)
	{
		const int l = (int)hgap_.size();
		for (int i = 0; i < l; ++i) {
			set_channel(hgap_[i], c, ::DISPATCH_ARCH::ScoreTraits<_sv>::zero_score());
			set_channel(score_[i], c, ::DISPATCH_ARCH::ScoreTraits<_sv>::zero_score());
		}
		set_channel(score_[l], c, ::DISPATCH_ARCH::ScoreTraits<_sv>::zero_score());
	}

	int cols() const {
		return cols_;
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
	int rows_, cols_;
};

#ifndef __APPLE__
template<typename _sv> thread_local MemBuffer<_sv> TracebackVectorMatrix<_sv>::hgap_;
template<typename _sv> thread_local MemBuffer<_sv> TracebackVectorMatrix<_sv>::score_;
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
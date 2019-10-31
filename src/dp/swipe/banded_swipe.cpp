/****
DIAMOND protein aligner
Copyright (C) 2013-2019 Benjamin Buchfink <buchfink@gmail.com>

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

#include <algorithm>
#include "../dp.h"
#include "swipe.h"
#include "target_iterator.h"
#include "../../util/data_structures/mem_buffer.h"

namespace DP { namespace BandedSwipe {
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
			*score_ptr_ = x;
		}
		inline void set_zero() {}
		_sv *hgap_ptr_, *score_ptr_;
	};
	Matrix(int band, size_t cols):
		band_(band)
	{
		hgap_.resize(band + 1);
		score_.resize(band);
		std::fill(hgap_.begin(), hgap_.end(), _sv());
		std::fill(score_.begin(), score_.end(), _sv());
		
	}
	inline ColumnIterator begin(int offset, int col)
	{
		return ColumnIterator(&hgap_[offset], &score_[offset]);
	}
	int band() const {
		return band_;
	}
private:
	int band_;
	static thread_local MemBuffer<_sv> hgap_, score_;
};

template<typename _sv>
struct TracebackMatrix
{

	typedef typename ScoreTraits<_sv>::Score Score;

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
			return *(score_ - band_ * ScoreTraits<_sv>::CHANNELS);
		}
		void walk_diagonal()
		{
			score_ -= band_ * ScoreTraits<_sv>::CHANNELS;
			--i;
			--j;
			assert(i >= -1 && j >= -1);
		}
		pair<Edit_operation, int> walk_gap(int d0, int d1)
		{
			const int i0 = std::max(d0 + j, 0), j0 = std::max(i - d1, -1);
			const Score *h = score_ - (band_ - 1) * ScoreTraits<_sv>::CHANNELS, *h0 = score_ - (j - j0) * (band_ - 1) * ScoreTraits<_sv>::CHANNELS;
			const Score *v = score_ - ScoreTraits<_sv>::CHANNELS, *v0 = score_ - (i - i0 + 1) * ScoreTraits<_sv>::CHANNELS;
			const Score score = this->score();
			const Score e = score_matrix.gap_extend();
			Score g = score_matrix.gap_open() + e;
			int l = 1;
			while (v > v0 && h > h0) {
				if (score + g == *h) {
					walk_hgap(h, l);
					return std::make_pair(op_deletion, l);
				}
				else if (score + g == *v) {
					walk_vgap(v, l);
					return std::make_pair(op_insertion, l);
				}
				h -= (band_ - 1) * ScoreTraits<_sv>::CHANNELS;
				v -= ScoreTraits<_sv>::CHANNELS;
				++l;
				g += e;
			}
			while (v > v0) {
				if (score + g == *v) {
					walk_vgap(v, l);
					return std::make_pair(op_insertion, l);
				}
				v -= ScoreTraits<_sv>::CHANNELS;
				++l;
				g += e;
			}
			while (h > h0) {
				if (score + g == *h) {
					walk_hgap(h, l);
					return std::make_pair(op_deletion, l);
				}
				h -= (band_ - 1) * ScoreTraits<_sv>::CHANNELS;
				++l;
				g += e;
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
		for (int i = i_; i < i1; ++i, s += ScoreTraits<_sv>::CHANNELS)
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

private:

	const size_t band_;
	static thread_local MemBuffer<_sv> hgap_, score_;

};

template<typename _sv> thread_local MemBuffer <_sv> Matrix<_sv>::hgap_;
template<typename _sv> thread_local MemBuffer<_sv> Matrix<_sv>::score_;
template<typename _sv> thread_local MemBuffer<_sv> TracebackMatrix<_sv>::hgap_;
template<typename _sv> thread_local MemBuffer<_sv> TracebackMatrix<_sv>::score_;

template<typename _sv, typename _traceback>
struct MatrixTag
{};

template<typename _sv>
struct MatrixTag<_sv, Traceback>
{
	typedef TracebackMatrix<_sv> Type;
};

template<typename _sv>
struct MatrixTag<_sv, ScoreOnly>
{
	typedef Matrix<_sv> Type;
};

template<typename _sv>
void traceback(const sequence &query, Frame frame, const TracebackMatrix<_sv> &dp, DpTarget &target, typename ScoreTraits<_sv>::Score max_score, int max_col, int channel, int i0, int i1)
{
	typedef typename ScoreTraits<_sv>::Score Score;
	const int j0 = i1 - (target.d_end - 1), d1 = target.d_end, d0 = target.d_begin;
	typename TracebackMatrix<_sv>::TracebackIterator it(dp.traceback(max_col + 1, i0 + max_col, j0 + max_col, (int)query.length(), channel, max_score));
	target.out->emplace_back();

	Hsp &out = target.out->back();
	out.score = target.score = ScoreTraits<_sv>::int_score(max_score);
	out.transcript.reserve(size_t(out.score * config.transcript_len_estimate));

	out.frame = frame.index();
	out.query_range.end_ = it.i + 1;
	out.subject_range.end_ = it.j + 1;
	
	while (it.score() > ScoreTraits<_sv>::zero_score()) {
		const Letter q = query[it.i], s = target.seq[it.j];
		const Score m = score_matrix(q, s), score = it.score();
		if (score == it.diag() + m) {
			out.push_match(q, s, m > (Score)0);
			it.walk_diagonal();
		} else {
			const pair<Edit_operation, int> g(it.walk_gap(d0, d1));
			out.push_gap(g.first, g.second, &target.seq[it.j + g.second]);
		}
	}

	out.query_range.begin_ = it.i + 1;
	out.subject_range.begin_ = it.j + 1;
	out.transcript.reverse();
	out.transcript.push_terminator();
}

template<typename _sv>
void traceback(const sequence &query, Frame frame, const Matrix<_sv> &dp, DpTarget &target, typename ScoreTraits<_sv>::Score max_score, int max_col, int channel, int i0, int i1)
{
	target.out->emplace_back();
	Hsp &out = target.out->back();
	out.score = target.score = ScoreTraits<_sv>::int_score(max_score);
	out.frame = frame.index();
}

template<typename _sv, typename _traceback>
void swipe(const sequence &query, Frame frame, vector<DpTarget>::iterator subject_begin, vector<DpTarget>::iterator subject_end)
{
	typedef typename ScoreTraits<_sv>::Score Score;
	typedef typename MatrixTag<_sv, _traceback>::Type Matrix;

	assert(subject_end - subject_begin <= ScoreTraits<_sv>::CHANNELS);
	const int qlen = (int)query.length();

	int band = 0;
	for (vector<DpTarget>::const_iterator j = subject_begin; j < subject_end; ++j)
		band = std::max(band, j->d_end - j->d_begin);

	int i1 = INT_MAX;
	for (vector<DpTarget>::iterator j = subject_begin; j < subject_end; ++j) {
		j->d_begin = j->d_end - band;
		i1 = std::min(i1, std::max(j->d_end - 1, 0));
	}
	int i0 = i1 + 1 - band;

	TargetIterator<ScoreTraits<_sv>::CHANNELS> targets(subject_begin, subject_end, i1, qlen);
	Matrix dp(band, targets.cols);

	const _sv open_penalty(static_cast<char>(score_matrix.gap_open() + score_matrix.gap_extend())),
		extend_penalty(static_cast<char>(score_matrix.gap_extend()));
	SwipeProfile<_sv> profile;

	Score best[ScoreTraits<_sv>::CHANNELS];
	int max_col[ScoreTraits<_sv>::CHANNELS];
	for (int i = 0; i < ScoreTraits<_sv>::CHANNELS; ++i) {
		best[i] = ScoreTraits<_sv>::zero_score();
		max_col[i] = 0;
	}

	int j = 0;
	while (targets.active.size() > 0) {
		const int i0_ = std::max(i0, 0), i1_ = std::min(i1, qlen - 1);
		if (i0_ > i1_)
			break;
		typename Matrix::ColumnIterator it(dp.begin(i0_ - i0, j));
		_sv vgap = _sv(), hgap = _sv(), col_best = _sv();
		if (i0_ - i0 > 0)
			it.set_zero();

		profile.set(targets.get());
		for (int i = i0_; i <= i1_; ++i) {
			hgap = it.hgap();
			const _sv next = cell_update<_sv>(it.diag(), profile.get(query[i]), extend_penalty, open_penalty, hgap, vgap, col_best);
			it.set_hgap(hgap);
			it.set_score(next);
			++it;
		}

		Score col_best_[ScoreTraits<_sv>::CHANNELS];
		store_sv(col_best, col_best_);
		for (int i = 0; i < targets.active.size();) {
			int channel = targets.active[i];
			if (!targets.inc(channel))
				targets.active.erase(i);
			else
				++i;
			if (col_best_[channel] > best[channel]) {
				best[channel] = col_best_[channel];
				max_col[channel] = j;
			}
		}
		++i0;
		++i1;
		++j;
	}

	for (int i = 0; i < targets.n_targets; ++i) {
		if (best[i] < ScoreTraits<_sv>::max_score()) {
			subject_begin[i].overflow = false;
			traceback<_sv>(query, frame, dp, subject_begin[i], best[i], max_col[i], i, i0 - j, i1 - j);
		}
		else
			subject_begin[i].overflow = true;
	}
}

template void swipe<score_vector<int16_t>, Traceback>(const sequence&, Frame, vector<DpTarget>::iterator, vector<DpTarget>::iterator);
template void swipe<score_vector<int16_t>, ScoreOnly>(const sequence&, Frame, vector<DpTarget>::iterator, vector<DpTarget>::iterator);

}}}
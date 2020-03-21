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
#include <utility>
#include <list>
#include <limits.h>
#include "../dp.h"
#include "swipe.h"
#include "target_iterator.h"
#include "../../util/data_structures/mem_buffer.h"
#include "../score_vector_int16.h"
#include "../../util/math/integer.h"
#include "../score_vector_int8.h"
#include "../../basic/config.h"
#include "../util/data_structures/range_partition.h"

using std::list;
using std::pair;

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
	static thread_local MemBuffer<_sv> hgap_, score_;
private:
	int band_;	
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
				h -= (band_ - 1) * ScoreTraits<_sv>::CHANNELS;
				v -= ScoreTraits<_sv>::CHANNELS;
				++l;
				score += e;
			}
			while (v > v0) {
				if (score == *v) {
					walk_vgap(v, l);
					return std::make_pair(op_insertion, l);
				}
				v -= ScoreTraits<_sv>::CHANNELS;
				++l;
				score += e;
			}
			while (h > h0) {
				if (score == *h) {
					walk_hgap(h, l);
					return std::make_pair(op_deletion, l);
				}
				h -= (band_ - 1) * ScoreTraits<_sv>::CHANNELS;
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

	static thread_local MemBuffer<_sv> hgap_, score_;

private:

	const size_t band_;
	

};

template<typename _sv> thread_local MemBuffer<_sv> Matrix<_sv>::hgap_;
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

template<typename _score>
_score add_cbs(_score x, int8_t b) {
	return x + _score(b);
}

template<typename _score>
_score add_cbs(_score x, void *b) {
	return x;
}

template<typename _sv, typename _cbs>
Hsp traceback(const sequence &query, Frame frame, _cbs bias_correction, const TracebackMatrix<_sv> &dp, const DpTarget &target, int d_begin, typename ScoreTraits<_sv>::Score max_score, int max_col, int channel, int i0, int i1)
{
	typedef typename ScoreTraits<_sv>::Score Score;
	const int j0 = i1 - (target.d_end - 1), d1 = target.d_end;
	typename TracebackMatrix<_sv>::TracebackIterator it(dp.traceback(max_col + 1, i0 + max_col, j0 + max_col, (int)query.length(), channel, max_score));
	
	Hsp out;
	out.swipe_target = target.target_idx;
	out.score = ScoreTraits<_sv>::int_score(max_score);
	out.transcript.reserve(size_t(out.score * config.transcript_len_estimate));

	out.frame = frame.index();
	out.query_range.end_ = it.i + 1;
	out.subject_range.end_ = it.j + 1;
	
	while (it.score() > ScoreTraits<_sv>::zero_score()) {
		const Letter q = query[it.i], s = target.seq[it.j];
		const Score m = score_matrix(q, s), score = it.score();
		const Score m2 = add_cbs(m, bias_correction[it.i]);
		if (score == saturated_add(it.diag(), m2)) {
			out.push_match(q, s, m > (Score)0);
			it.walk_diagonal();
		} else {
			const pair<Edit_operation, int> g(it.walk_gap(d_begin, d1));
			out.push_gap(g.first, g.second, &target.seq[it.j + g.second]);
		}
	}

	out.query_range.begin_ = it.i + 1;
	out.subject_range.begin_ = it.j + 1;
	out.transcript.reverse();
	out.transcript.push_terminator();
	return out;
}

template<typename _sv, typename _cbs>
Hsp traceback(const sequence &query, Frame frame, _cbs bias_correction, const Matrix<_sv> &dp, const DpTarget &target, int d_begin, typename ScoreTraits<_sv>::Score max_score, int max_col, int channel, int i0, int i1)
{
	Hsp out;
	out.swipe_target = target.target_idx;
	out.score = ScoreTraits<_sv>::int_score(max_score);
	out.frame = frame.index();
	out.query_range = interval(target.d_begin, target.d_end);
	out.subject_range = interval(target.j_begin, target.j_end);
	return out;
}

template<typename _traceback>
bool realign(const Hsp &hsp, const DpTarget &dp_target) {
	return false;
}

template<>
bool realign<Traceback>(const Hsp &hsp, const DpTarget &dp_target) {
	return hsp.subject_range.begin_ - config.min_realign_overhang > dp_target.j_begin || hsp.subject_range.end_ + config.min_realign_overhang  < dp_target.j_end;
}

template<typename _sv, typename _traceback, typename _cbs>
list<Hsp> swipe(
	const sequence &query,
	Frame frame,
	vector<DpTarget>::const_iterator subject_begin,
	vector<DpTarget>::const_iterator subject_end,
	_cbs composition_bias,
	int score_cutoff,
	vector<DpTarget> &overflow,
	Statistics &stat)
{
	typedef typename ScoreTraits<_sv>::Score Score;
	typedef typename MatrixTag<_sv, _traceback>::Type Matrix;
 	constexpr int CHANNELS = ScoreTraits<_sv>::CHANNELS;

	assert(subject_end - subject_begin <= CHANNELS);
	const int qlen = (int)query.length();

	int band = 0;
	for (vector<DpTarget>::const_iterator j = subject_begin; j < subject_end; ++j)
		band = std::max(band, j->d_end - j->d_begin);

	int i1 = INT_MAX, d_begin[CHANNELS];
	const int target_count = int(subject_end - subject_begin);
#ifdef STRICT_BAND
	int band_offset[CHANNELS];
#endif
	for (int i = 0; i < target_count; ++i) {
		d_begin[i] = subject_begin[i].d_end - band;
#ifdef STRICT_BAND
		band_offset[i] = subject_begin[i].d_begin - d_begin[i];
#endif		
		i1 = std::min(i1, std::max(subject_begin[i].d_end - 1, 0));
	}
	int i0 = i1 + 1 - band;
#ifdef STRICT_BAND
	RangePartition<CHANNELS, Score> band_parts(band_offset, target_count, band);
#endif
	
	TargetIterator<CHANNELS> targets(subject_begin, subject_end, i1, qlen, d_begin);
	Matrix dp(band, targets.cols);

	const _sv open_penalty(static_cast<char>(score_matrix.gap_open() + score_matrix.gap_extend())),
		extend_penalty(static_cast<char>(score_matrix.gap_extend()));
	SwipeProfile<_sv> profile;

	Score best[CHANNELS];
	int max_col[CHANNELS];
	std::fill(best, best + CHANNELS, ScoreTraits<_sv>::zero_score());
	std::fill(max_col, max_col + CHANNELS, 0);
	
	int j = 0;
	while (targets.active.size() > 0) {
		const int i0_ = std::max(i0, 0), i1_ = std::min(i1, qlen - 1) + 1;
		if (i0_ >= i1_)
			break;
		typename Matrix::ColumnIterator it(dp.begin(i0_ - i0, j));
		_sv vgap = _sv(), hgap = _sv(), col_best = _sv();
		if (i0_ - i0 > 0)
			it.set_zero();

		profile.set(targets.get(Score()));

#ifdef STRICT_BAND
		for (int part = 0; part < band_parts.count(); ++part) {
			const int i_begin = std::max(i0 + band_parts.begin(part), i0_);
			const int i_end = std::min(i0 + band_parts.end(part), i1_);
			const _sv target_mask = load_sv(band_parts.mask(part));
			for (int i = i_begin; i < i_end; ++i) {
#else
			for (int i = i0_; i < i1_; ++i) {
#endif
				hgap = it.hgap();
				_sv next;
				_sv match_scores = profile.get(query[i]);
#ifdef STRICT_BAND
				match_scores += target_mask;
#endif
				next = swipe_cell_update<_sv>(it.diag(), match_scores, composition_bias[i], extend_penalty, open_penalty, hgap, vgap, col_best);

				it.set_hgap(hgap);
				it.set_score(next);
				++it;
			}
#ifdef STRICT_BAND
		}
#endif

		Score col_best_[CHANNELS];
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

	list<Hsp> out;
	int realign = 0;
	for (int i = 0; i < targets.n_targets; ++i) {
		if (best[i] < ScoreTraits<_sv>::max_score()) {
			if (ScoreTraits<_sv>::int_score(best[i]) >= score_cutoff) {
				out.push_back(traceback<_sv>(query, frame, composition_bias, dp, subject_begin[i], d_begin[i], best[i], max_col[i], i, i0 - j, i1 - j));
				if ((config.max_hsps == 0 || config.max_hsps > 1) && !config.no_swipe_realign
					&& ::DP::BandedSwipe::DISPATCH_ARCH::realign<_traceback>(out.back(), subject_begin[i]))
					realign |= 1 << i;
			}
		}
		else
			overflow.push_back(subject_begin[i]);
	}
	if (realign) {
		stat.inc(Statistics::SWIPE_REALIGN);
		vector<vector<Letter>> seqs;
		vector<DpTarget> realign_targets;
		for (int i = 0; i < targets.n_targets; ++i) {
			if ((realign & (1 << i)) == 0)
				continue;
			seqs.push_back(subject_begin[i].seq.copy());
			realign_targets.push_back(subject_begin[i]);
			realign_targets.back().seq = sequence(seqs.back());
			for (const Hsp& hsp : out)
				if (hsp.swipe_target == subject_begin[i].target_idx)
					realign_targets.back().seq.mask(hsp.subject_range);
		}
		vector<DpTarget> overflow;
		out.splice(out.end(), swipe<_sv, _traceback>(query, frame, realign_targets.begin(), realign_targets.end(), composition_bias, score_cutoff, overflow, stat));
	}
	return out;
}

#ifdef __SSE4_1__
template list<Hsp> swipe<score_vector<int8_t>, Traceback, const int8_t*>(const sequence&, Frame, vector<DpTarget>::const_iterator, vector<DpTarget>::const_iterator, const int8_t*, int, vector<DpTarget>&, Statistics&);
template list<Hsp> swipe<score_vector<int8_t>, ScoreOnly, const int8_t*>(const sequence&, Frame, vector<DpTarget>::const_iterator, vector<DpTarget>::const_iterator, const int8_t*, int, vector<DpTarget>&, Statistics&);
#endif
#ifdef __SSE2__
template list<Hsp> swipe<score_vector<int16_t>, Traceback, const int8_t*>(const sequence&, Frame, vector<DpTarget>::const_iterator, vector<DpTarget>::const_iterator, const int8_t*, int, vector<DpTarget>&, Statistics&);
template list<Hsp> swipe<score_vector<int16_t>, ScoreOnly, const int8_t*>(const sequence&, Frame, vector<DpTarget>::const_iterator, vector<DpTarget>::const_iterator, const int8_t*, int, vector<DpTarget>&, Statistics&);
#endif
template list<Hsp> swipe<int32_t, Traceback, const int8_t*>(const sequence&, Frame, vector<DpTarget>::const_iterator, vector<DpTarget>::const_iterator, const int8_t*, int, vector<DpTarget>&, Statistics&);
template list<Hsp> swipe<int32_t, ScoreOnly, const int8_t*>(const sequence&, Frame, vector<DpTarget>::const_iterator, vector<DpTarget>::const_iterator, const int8_t*, int, vector<DpTarget>&, Statistics&);

#ifdef __SSE4_1__
template list<Hsp> swipe<score_vector<int8_t>, Traceback, NoCBS>(const sequence&, Frame, vector<DpTarget>::const_iterator, vector<DpTarget>::const_iterator, NoCBS, int, vector<DpTarget>&, Statistics&);
template list<Hsp> swipe<score_vector<int8_t>, ScoreOnly, NoCBS>(const sequence&, Frame, vector<DpTarget>::const_iterator, vector<DpTarget>::const_iterator, NoCBS, int, vector<DpTarget>&, Statistics&);
#endif
#ifdef __SSE2__
template list<Hsp> swipe<score_vector<int16_t>, Traceback, NoCBS>(const sequence&, Frame, vector<DpTarget>::const_iterator, vector<DpTarget>::const_iterator, NoCBS, int, vector<DpTarget>&, Statistics&);
template list<Hsp> swipe<score_vector<int16_t>, ScoreOnly, NoCBS>(const sequence&, Frame, vector<DpTarget>::const_iterator, vector<DpTarget>::const_iterator, NoCBS, int, vector<DpTarget>&, Statistics&);
#endif
template list<Hsp> swipe<int32_t, Traceback, NoCBS>(const sequence&, Frame, vector<DpTarget>::const_iterator, vector<DpTarget>::const_iterator, NoCBS, int, vector<DpTarget>&, Statistics&);
template list<Hsp> swipe<int32_t, ScoreOnly, NoCBS>(const sequence&, Frame, vector<DpTarget>::const_iterator, vector<DpTarget>::const_iterator, NoCBS, int, vector<DpTarget>&, Statistics&);

}}}
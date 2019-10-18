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
		_sv *hgap_ptr_, *score_ptr_;
	};
	Matrix(int band):
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
	struct ColumnIterator
	{
		ColumnIterator(_sv* hgap_front, _sv* score_front, _sv* hgap_front1, _sv* score_front1) :
			hgap_ptr_(hgap_front),
			score_ptr_(score_front),
			hgap_ptr1_(hgap_front1),
			score_ptr1_(score_front1)
		{ }
		inline void operator++()
		{
			++hgap_ptr_; ++score_ptr_; ++hgap_ptr1_; ++score_ptr1_;
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
			*hgap_ptr1_ = x;
		}
		inline void set_score(const _sv& x)
		{
			*score_ptr1_ = x;
		}
		void set_zero()
		{
			(score_ptr1_ - 1)->zero();
		}
		_sv *hgap_ptr_, *score_ptr_, *hgap_ptr1_, *score_ptr1_;
	};
	TracebackMatrix(size_t band, size_t cols) :
		band_(band)
	{
		hgap_.resize((band + 1) * (cols + 1));
		score_.resize(band * (cols + 1));
		size_t i = 0;
		_sv z = _sv()
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
	static thread_local MemBuffer<_sv> hgap_, score_;
};

template<typename _sv> thread_local MemBuffer <_sv> Matrix<_sv>::hgap_;
template<typename _sv> thread_local MemBuffer<_sv> Matrix<_sv>::score_;
template<typename _sv> thread_local MemBuffer<_sv> TracebackMatrix<_sv>::hgap_;
template<typename _sv> thread_local MemBuffer<_sv> TracebackMatrix<_sv>::score_;

template<typename _sv>
void traceback(const sequence &query, Strand strand, int dna_len, const Matrix<_sv> &dp, DpTarget &target, typename ScoreTraits<_sv>::Score max_score, int max_col, int channel, int i0, int i1, bool parallel)
{
	if (parallel)
		target.tmp = new Hsp();
	else
		target.out->emplace_back();
	Hsp &out = parallel ? *target.tmp : target.out->back();

	const int j0 = i1 - (target.d_end - 1);
	out.score = target.score = ScoreTraits<_sv>::int_score(max_score);
	out.query_range.end_ = std::min(i0 + max_col + (int)dp.band() / 3 / 2, (int)query.length());
	out.query_range.begin_ = std::max(out.query_range.end_ - (j0 + max_col), 0);
	out.frame = strand == FORWARD ? 0 : 3;
	out.query_source_range = TranslatedPosition::absolute_interval(TranslatedPosition(out.query_range.begin_, Frame(out.frame)), TranslatedPosition(out.query_range.end_, Frame(out.frame)), dna_len);
}

template<typename _sv>
void swipe(const sequence &query, vector<DpTarget>::iterator subject_begin, vector<DpTarget>::iterator subject_end)
{
	typedef typename ScoreTraits<_sv>::Score Score;

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
	Matrix<_sv> dp(band);

	const _sv open_penalty(static_cast<char>(score_matrix.gap_open() + score_matrix.gap_extend())),
		extend_penalty(static_cast<char>(score_matrix.gap_extend()));
	_sv best = _sv();
	SwipeProfile<_sv> profile;

	int j = 0;
	while (targets.active.size() > 0) {
		const int i0_ = std::max(i0, 0), i1_ = std::min(i1, qlen - 1);
		if (i0_ > i1_)
			break;
		typename Matrix<_sv>::ColumnIterator it(dp.begin(i0_ - i0, j));
		_sv vgap = _sv(), hgap = _sv();

		profile.set(targets.get());
		for (int i = i0_; i <= i1_; ++i) {
			hgap = it.hgap();
			const _sv next = cell_update<_sv>(it.diag(), profile.get(query[i]), extend_penalty, open_penalty, hgap, vgap, best);
			it.set_hgap(hgap);
			it.set_score(next);
			++it;
		}

		for (int i = 0; i < targets.active.size();) {
			int channel = targets.active[i];
			if (!targets.inc(channel))
				targets.active.erase(i);
			else
				++i;
		}
		++i0;
		++i1;
		++j;
	}

	Score max_score[ScoreTraits<_sv>::CHANNELS];
	best.store(max_score);
	for (int i = 0; i < targets.n_targets; ++i) {
		subject_begin[i].overflow = false;
		traceback<_sv>(query, FORWARD, (int)query.length(), dp, subject_begin[i], max_score[i], 0, i, i0 - j, i1 - j, false);
	}
}

template<typename _sv>
void swipe_targets(const sequence &query,
	vector<DpTarget>::iterator begin,
	vector<DpTarget>::iterator end)
{
	for (vector<DpTarget>::iterator i = begin; i < end; i += ScoreTraits<_sv>::CHANNELS) {
		/*if (!overflow_only || i->overflow) {
			if (score_only || config.disable_traceback)
				banded_3frame_swipe<_sv, ScoreOnly>(query, strand, i, i + std::min(vector<DpTarget>::iterator::difference_type(ScoreTraits<_sv>::CHANNELS), end - i), stat, parallel);
			else
				banded_3frame_swipe<_sv, Traceback>(query, strand, i, i + std::min(vector<DpTarget>::iterator::difference_type(ScoreTraits<_sv>::CHANNELS), end - i), stat, parallel);
		}*/
		swipe<_sv>(query, i, i + std::min(vector<DpTarget>::iterator::difference_type(ScoreTraits<_sv>::CHANNELS), end - i));
	}
}

void swipe(const sequence &query, vector<DpTarget>::iterator target_begin, vector<DpTarget>::iterator target_end)
{
#ifdef __SSE2__
	std::stable_sort(target_begin, target_end);
	swipe_targets<score_vector<int16_t>>(query, target_begin, target_end);
#endif
}

}}